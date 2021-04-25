library(tidyverse)
library(data.table)
library(DESeq2)
library(pacman)

theme_set(theme_classic())

setwd("~/bone_data_analysis")

# set up the data and metadata
bone_data = fread("htseq_all_counts.txt") %>% as.data.frame()
rownames(bone_data) = bone_data$Feature
bone_data = bone_data[,-1]

# bone_data = fread("alldata_raw_counts.txt")[, -c(1:8)] %>% as.data.frame()
colData = fread("bone_coldata.txt", stringsAsFactors = T) %>% as.data.frame()
rownames(colData) = colData$`Sample name`
bone_data = bone_data[, rownames(colData)]

# check
all(rownames(colData) == colnames(bone_data))

trap_cre_mask = colData$Cre == "Trap"
trap_cre_mask[is.na(trap_cre_mask)] = F
bone_data_trap = bone_data[, trap_cre_mask]
colData_trap = colData[trap_cre_mask,]

# check
all(rownames(colData_trap) == colnames(bone_data_trap))

# colData_trap[] = lapply(colData_trap, factor)
# resected_cat = as.factor(ifelse(colData$Dpi == "1 week", 1, 0))
#bone_data_trap = bone_data_trap[, -c(1:4)]

dds = DESeqDataSetFromMatrix(bone_data_trap, colData_trap, design = ~Dpi)
DE_obj = DESeq(dds)
res = results(DE_obj)
resLFC = lfcShrink(DE_obj, coef = "Dpi_unresected_vs_1.week", type = "apeglm")
MA_plot = plotMA(resLFC)

# PCA 
vsd <- vst(DE_obj, blind=FALSE)
PCA_plot = plotPCA(vsd, intgroup = "Dpi")

sig_thresh = 0.05
volc_plot_data = data.frame(res)
volc_plot_data$signif = ifelse(!is.na(volc_plot_data$padj) && volc_plot_data$padj <= sig_thresh, "sig", "nsig")

cols = c("red", "black")
names(cols) = c("sig", "nsig")
volc_plot = ggplot(data = volc_plot_data, aes(x = log2FoldChange, y = -log10(padj), color = signif)) + 
  geom_point() + 
  geom_hline(yintercept = -log10(sig_thresh), col = "red") +
  scale_color_manual(values = cols)

print(cat("there are ", sum(volc_plot_data$signif == "sig"), " DEGs"))

p_load(pheatmap)
top50genes = rownames(volc_plot_data)[order(volc_plot_data$padj)[1:50]]
interest_genes = c(top50genes, "Flt3", "Cx3cr1", "Csf1r")
heatmap_data = scale(t(bone_data_trap[which(rownames(bone_data_trap) %in% interest_genes), ]))

annotation_col = data.frame(
  Dpr = factor(colData_trap$Dpi)
)
rownames(annotation_col) = rownames(heatmap_data)

pheatmap(heatmap_data, annotation_row = annotation_col)
