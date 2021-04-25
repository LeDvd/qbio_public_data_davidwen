library(pacman)

library(maftools)
library(TCGAbiolinks)
library(S4Vectors)
p_load(tidyverse)

setwd("~/qbio_public_data_davidwen/")

barcodes_rnaseq <- c("TCGA-BH-A0DG-01A-21R-A12P-07","TCGA-A2-A0YF-01A-21R-A109-07",
         "TCGA-AN-A04A-01A-21R-A034-07","TCGA-AR-A1AP-01A-11R-A12P-07",
             "TCGA-A2-A0T3-01A-21R-A115-07", "TCGA-E2-A154-01A-11R-A115-07")

# mutation <- GDCquery_Maf(tumor = "BRCA", save.csv=TRUE, pipeline="mutect2")
maf_df = read.csv("GDCdata/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.csv")
maf = read.maf(maf_df)

# # plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
# oncoplot(maf, top=10)
# somaticInteractions(maf, top = 25, pvalue = c(0.05, 0.1))
# 
# odrive = oncodrive(maf, minMut = 5, pvalMethod = "zscore")
# plotOncodrive(res = odrive, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
# 
# prog_geneset = survGroup(maf = maf, top = 20, geneSetSize = 2, time = "days_to_last_followup", Status = "Overall_Survival_Status", verbose = FALSE)
# 
# library("mclust")
# het = inferHeterogeneity(maf, tsb = 'TCGA-AN-A046-01A-21W-A050-09')
# plotClusters(het)

p_load(vegan)
p_load(cluster)

mtx = mutCountMatrix(maf, includeSyn = F, countOnly = NULL, removeNonMutated = F)
mtx_b = t(mtx) %>%
    apply(2, function(x) ifelse(x > 0, 1, x))
gower <- daisy(mtx_b, metric = c("gower"))
gower.aggl.clust<- hclust (gower, method = "complete")
# plot(gower.aggl.clust, cex = 0.6, main = "Agglomerative, complete linkages")

# hcd = as.dendrogram(gower.aggl.clust)
# plot(hcd, type="rectangle", ylab="Height")

p_load(ggdendro)

# ggdendrogram(gower.aggl.clust, labels=F, sub=F, rotate=T)

p_load(ConsensusClusterPlus)
results = ConsensusClusterPlus(gower, maxK=6, reps=50, pItem=0.8)
icl = calcICL(results)
