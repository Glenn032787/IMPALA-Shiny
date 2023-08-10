# Long POG IMPALA dashboard

This is an RShiny dashboard that displays allele specific expression (ASE) data for the Long POG cohort. The dashboard shows a table of the ASE data of all phased genes in each sample in the Long POG cohort which can be filtered down based on major expressing allele frequency and various biological mechanisms. Various summary figures is generated based on the filtered data as well as the option to download the data. 

### Data avalible
The data can be found here:

### Methods
Using the IMPALA pipeline, we obtained ASE data from 174 cancer samples. Oxford Nanopore long-reads data was available for this cohort, which was used to obtain haplotype information through genomic phasing. Combining RNA-seq data with haplotype information of SNVs, IMPALA pipeline was able to call ASE data with a higher accuracy. Additionally, ASE data is intersected with CNV data from ploidetect, allelic methylation data from NanoMethPhase and somatic data to explain the mechanism for the ASE. 

### Accessing it online
This dashboard be found on this website.


### Running manually
This app can be run manually by cloning this repository and using `runApp` command in R terminal.The 

```
git clone git@github.com:Glenn032787/IMPALA-Shiny.git
cd IMPALA-Shiny

R # Start R terminal 
> library(shiny)
> shiny::runApp("app.R")
Listening on http://127.0.0.1:4305
```
