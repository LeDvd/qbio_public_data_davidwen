library(SummarizedExperiment) #look up how a SummarizedExperiment is built
library(TCGAbiolinks)

setwd("/Users/Turnip/qbio_public_data_davidwen")

clinical_file <- read.csv("data/tcga_brca_six_example_clinical.csv")
barcodes <- as.character( clinical_file$barcode )
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  barcode = c(barcodes))
GDCdownload(query)
data <- GDCprepare(query)
str(data) #use this line if you are typing in the command line
