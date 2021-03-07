library(pacman)
p_load(TCGAbiolinks)
p_load(DESeq2)
p_load(tidyverse)

# REFERENCE: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# setwd("~/qbio_public_data_davidwen/R_code")

# makes the plots pretty
theme_set(theme_bw())

# clinical_file <- read.csv("../data/tcga_brca_six_example_clinical.csv")
# barcodes <- as.character( clinical_file$barcode )

# you know the drill
query = GDCquery(project = "TCGA-BRCA", 
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification",
                 workflow.type = "HTSeq - Counts")

# GDCdownload(query)
data = GDCprepare(query)

# need to remove NA values
data = data[, !is.na(data$age_at_index)]
data$age_cat = factor(ifelse(data$age_at_index >= 60, "old", 
                        ifelse(data$age_at_index <= 40, "young", "mid")),
                      levels = c("young", "mid", "old"))


# standard DESeq pipeline
dds = DESeqDataSet(data, design = ~ age_cat)
DE_obj = DESeq(dds)
saveRDS(DE_obj, file = "deseq.rds") # i dont wanna run this again

res = results(DE_obj)
resLFC = lfcShrink(DE_obj, coef = "age_cat_old_vs_young", type = "apeglm")
resOrdered = res[order(res$padj), ]

# check for biases
MA_plot = plotMA(resLFC)

# PCA with normalized results
vsd <- vst(DE_obj, blind=FALSE)
PCA_plot = plotPCA(vsd, intgroup = "age_cat")

# create a volcano plot
volc_plot_data = data.frame(res)
volc_plot_data$signif = ifelse(volc_plot_data$padj <= 0.05, "sig", "nsig")
cols = c("red", "black")
names(cols) = c("sig", "nsig")
volc_plot = ggplot(data = volc_plot_data, aes(x = log2FoldChange, y = -log10(padj), color = signif)) + 
               geom_point() + 
               geom_hline(yintercept = -log10(0.05), col = "red") +
               scale_color_manual(values = cols)

# save outputs
ggsave("/diff_analysis/PCA_plot.png", plot = PCA_plot, device = "png", width = 11, height = 8.5, units = "in")
ggsave("/diff_analysis/volc_plot.png", plot = volc_plot, device = "png", width = 11, height = 8.5, units = "in")
ggsave("/diff_analysis/MA_plot.png", plot = MA_plot, device = "png", width = 11, height = 8.5, units = "in")
