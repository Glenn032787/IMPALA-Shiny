library(shiny)
library(tidyverse)
library(DT)
library(shinyWidgets)

options(dplyr.summarise.inform = FALSE)
options(digits=2)

# Get ASE data into shiny app
aseData <- readRDS("/projects/glchang_prj/finalPOGdata/allASEdata_anonymize.RDS") %>%
    dplyr::mutate(methylation = case_when(
        methyl_state < 0 ~ "allele2Methyl",
        TRUE ~ "allele1Methyl"
    )) %>%
    dplyr::mutate(expression = case_when(
        allele1IsMajor ~ "alelle1Expression",
        TRUE ~ "alelle2Expression"
    )) %>%
    dplyr::mutate(aseResult = case_when(
        majorAlleleFrequency < 0.65 ~ "BAE",
        majorAlleleFrequency >= 0.65 ~ "ASE"
    )) %>%
    dplyr::mutate(majorAlleleFrequency = format(round(majorAlleleFrequency, 3), nsmall = 3)) %>%
    dplyr::mutate(padj = format(round(padj, 3), nsmall = 3))
    
pogSample <- c("All", unique(aseData$sample))
defaultColumn <- c("gene", "majorAlleleFrequency", 
                   "padj", "aseResults", "sample")
column <- setdiff(colnames(aseData), c(defaultColumn, "methylation", "expression", "aseResult"))
geneOfInterest <- readRDS("src/geneOfInterest.RDS")


# Define UI for application that draws a histogram
ui <- fluidPage(
   
    # Application title
    div(
        style={'padding: 15px'},
        img(src = "impala_logo.svg", align = "right", width = 150, height = 150),
        titlePanel("ASE data for Long POG cohort"),
        
        h5("This shiny app is used to explore Allele Specific Expression (ASE) data 
       from the Long POG cohort. Using Oxford Nanopore long reads for genome phasing 
       and Illumina RNA-seq data, the ", a("IMPALA pipeline", href="https://github.com/bcgsc/IMPALA/tree/master"),
           " was used to call ASE data in 174 cancer samples. Genes are annotated with 
       allelic CNV data, allelic methylation and nonsense mutation to explain the mechanism of ASE"),
        hr(),
        br(),
    
    
    fluidRow(
        sidebarPanel(
            width = 3,
            # Show number of ASE genes
            h4(textOutput("aseNum")),
            hr(),
            # Select Column
            pickerInput(
                "columnSelect",
                "Select columns to show",
                choices = column,
                multiple = TRUE
            ),
        
            # Filter for sample
            selectInput(
                "sample",
                label = "POG sample",
                choices = pogSample
            ),
            
            # Filter for gene
            textInput(
                "geneText",
                label = "Search for gene",
                value = ""
            ),
            
            # Filter for major allele frequency
            sliderInput("mafSlider",
                        "Major Allele Frequency",
                        min = 0.5, 
                        max = 1.0,
                        value = c(0.5, 1.0)),
            
            # Filter for allele
            radioButtons("allele1majorRadio",
                         "Filter major expressing haplotype",
                         choices = c("Both",
                                     "Haplotype 1", 
                                     "Haplotype 2"),
                         selected = "Both"),
            
            # Filter for gene type
            pickerInput("GeneOfInterest", 
                       "Special types of genes", 
                       choices = list("Chromosome X genes" = "chrX", 
                                      "Imprinting genes" = "imprinting", 
                                      "Olfactory genes" = "olfactory",
                                      "HLA genes" = "HLA", 
                                      "Tumor Suppressors" = "TumorSuppressor", 
                                      "Oncogenes" = "Oncogene"),
                        multiple = T),
            
            # Filter for CNV
            checkboxGroupInput("cnvCheckList", 
                               "Allelic CNV states", 
                               choices = list("balance" = "balance", 
                                              "imbalance" = "imbalance",
                                              "LOH" = "LOH",
                                              "No CNV data" = "NA"),
                               selected = c("balance", 
                                            "imbalance",
                                            "LOH",
                                            "NA")),
            
            # Filter for Methylation
            checkboxGroupInput("allelicMethyl", 
                               "Allelic Methylation", 
                               choices = list("Haplotype 1", 
                                              "Haplotype 2",
                                              "No allelic methylation"),
                               selected = c("Haplotype 1", 
                                            "Haplotype 2",
                                            "No allelic methylation")),
            
            # Filter for nonsense mutation
            checkboxGroupInput("nonsenseChecklist", 
                               "Nonsense mutation", 
                               choices = list("Haplotype 1" = 1, 
                                              "Haplotype 2" = 2,
                                              "No nonsense mutation" = "NA"),
                               selected = c(1,2,"NA")),
            
            # Filter for nonsense mutation
            checkboxGroupInput("somaticChecklist", 
                               "Somatic Mutation", 
                               choices = list("SNV", 
                                              "Indel")),
            
            # Filter for significance
            checkboxGroupInput("significance", 
                               "Filter for significant genes", 
                               choices = list("Only Significant")),
            hr(),
            downloadButton('downloadData', 'Download current table as csv')
        ),
        # Show a plot of the generated distribution
        mainPanel(
            
            tabsetPanel(type = "tabs",
                        tabPanel("Table",
                                 h4("Long POG data"),
                                 dataTableOutput("impala")),
                        tabPanel("Figures", 
                                 h4("Summary Figures"),
                                 br(),
                                 fluidRow(
                                     column(6, plotOutput("density")),
                                     column(6, plotOutput("general"))
                                 ),
                                 fluidRow(
                                     column(6, plotOutput("cnvBar")),
                                     column(6, plotOutput("methyl"))
                                 ))
                        )
        )
    )
)
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    methylFilter <- function(table, methyl){
        result <- tibble()
        if ("Haplotype 1" %in% methyl){
            result <- table %>%
                dplyr::filter(methyl_state >= 0) %>%
                bind_rows(result)
        } 
        
        if ("Haplotype 2" %in% methyl){
            result <- table %>%
                dplyr::filter(methyl_state < 0) %>%
                bind_rows(result)
        }
        
        if ("No allelic methylation" %in% methyl){
            result <- table %>%
                dplyr::filter(is.na(methyl_state)) %>%
                bind_rows(result)
        } 
        return(result)
    }

    
    table <- reactive({
        cnvCheckList <- replace(input$cnvCheckList, input$cnvCheckList=="NA", NA)
        
        tbl <- aseData %>%
            dplyr::filter(cnv_state %in% cnvCheckList) %>%
            dplyr::filter(grepl(input$geneText,gene) ) %>%
            dplyr::filter(majorAlleleFrequency >= input$mafSlider[1]) %>%
            dplyr::filter(majorAlleleFrequency <= input$mafSlider[2]) %>%
            dplyr::filter(stop_variant_allele %in% suppressWarnings(
                as.numeric(input$nonsenseChecklist))) 
        
        # Filter for Significance
        if (!is.null(input$significance)) {
            tbl <- tbl %>%
                dplyr::filter(padj <= 0.05)
        }
        
        # Filter for sample
        if (input$sample != "All") {
            tbl <- tbl %>%
                dplyr::filter(sample == input$sample)
        }
        
        # Filter for major expressing allele
        if (input$allele1majorRadio == "Haplotype 1") {
            tbl <- tbl %>%
                dplyr::filter(allele1IsMajor)
        } else if (input$allele1majorRadio == "Haplotype 2") {
            tbl <- tbl %>%
                dplyr::filter(!allele1IsMajor)
        }
        
        if (length(input$GeneOfInterest) != 0 ){
            specialGenes <- c()
            for (type in input$GeneOfInterest){
                specialGenes <- c(specialGenes, geneOfInterest[[type]])
            }
            tbl <- tbl %>%
                dplyr::filter(gene %in% specialGenes)
        }
        
        # Filter based on allelic methylation
        tbl <- methylFilter(tbl, input$allelicMethyl)
        
        # Filter somatic mutation
        if ("SNV" %in% input$somaticChecklist & "Indel" %in% input$somaticChecklist){
            tbl <- tbl %>%
                dplyr::filter(somaticSNV | somaticIndel)
        } else if ("SNV" %in% input$somaticChecklist) {
            tbl <- tbl %>%
                dplyr::filter(somaticSNV)
        }  else if ("Indel" %in% input$somaticChecklist) {
            tbl <- tbl %>%
                dplyr::filter(somaticIndel)
        } 
        
        tbl
    })
    
    output$impala <- DT::renderDataTable({
        table()  %>%
            dplyr::select(c(defaultColumn, input$columnSelect))
    })
    
    
    output$aseNum <- renderText(paste0(table() %>%
                                           dplyr::filter(aseResults == "ASE") %>%
                                           nrow() %>%
                                           format(nsmall=1, big.mark=","), 
                                       " ASE genes found"))
    output$cnvBar <- renderPlot(table() %>%
                                    dplyr::filter(padj < 0.05) %>%
                                    dplyr::filter(!is.na(cnv_state)) %>%
                                    ggplot(aes(aseResults)) + 
                                    geom_bar() + 
                                    facet_grid(~cnv_state) +
                                    ggtitle("Number of ASE and BAE gene in each CNV state")
    )
    
    output$general <- renderPlot(table() %>%
                                     group_by(sample, aseResult) %>%
                                     summarize(n=n()) %>%
                                     dplyr::mutate(prop = n/sum(n)) %>%
                                     ggplot(aes(aseResult, prop)) + 
                                     geom_boxplot() + 
                                     ggtitle("Proportion of ASE and BAE gene in each sample", 
                                             subtitle = "Major Allele Frequency above 0.65 = ASE") + 
                                     xlab("aseResults")
    )
    
    output$density <- renderPlot(table() %>%
                                     dplyr::mutate(majorAlleleFrequency = as.numeric(majorAlleleFrequency)) %>%
                                     ggplot(aes(majorAlleleFrequency)) + 
                                     geom_density() +
                                     geom_vline(xintercept = 0.65) +
                                     ggtitle("Density of Major Allele Frequency")
                                 )
    
    output$methyl <- renderPlot(table() %>%
                                dplyr::filter(!is.na(methyl_state)) %>%
                                dplyr::select(methylation, expression) %>%
                                group_by(methylation, expression) %>%
                                summarize(n=n()) %>%
                                ggplot(mapping = aes(x = methylation, y = expression)) +
                                geom_tile(aes(fill = n), colour = "white") +
                                geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, colour = "white") +
                                ggtitle("Allelic Methylation and Expression Contigency Table") +
                                theme_bw() + 
                                theme(legend.position = "none", 
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), 
                                      panel.border = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.ticks.y = element_blank(), 
                                      axis.title.x=element_blank(), 
                                      axis.title.y=element_blank())  
                                )
    output$downloadData <- downloadHandler(
        filename = function() {
            paste('IMPALA-data-', Sys.Date(), '.csv', sep='')
        },
        content = function(file) {
            write.csv(table(), file)
        },
        contentType = 'text/csv'
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
