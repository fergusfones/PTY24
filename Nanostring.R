setwd("/Users/fergusfones/Desktop/Nanostring/")
getwd()

library(standR)
library(SpatialExperiment)
library(ggplot2)
library(ggalluvial)
library(tidyverse)
library(readr)
library(ggpubr)
library(edgeR)
library(limma)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GeomxTools")
library(GeomxTools)

metaData <- read.csv("preProcessedMetaData_filtered.csv",header = T)
countData <- read.csv("preProcessedCountData_filtered.csv",header = T,)
pathology <- read.csv("PathologyPG5.csv", header = T,)
rownames(pathology) <- pathology$Histology.no.
rownames(metaData) <- metaData$SegmentDisplayName
colnames(countData)[14:108] <- gsub("X","",colnames(countData)[14:108])

duplicates <- c("D830030K20Rik", "Gm10406", "LOC118568634")
subCountData <- countData[-which(countData$TargetName %in% duplicates),]  
rownames(subCountData) <- subCountData$ProbeDisplayName

spe <- readGeoMx(countFile = subCountData[,c(4,14:108)],sampleAnnoFile = metaData,
                 featureAnnoFile = subCountData[,-c(14:108)],rmNegProbe = T) 

# general:
head(metaData)[, 1:5]
head(countData)[,1:5]
head(pathology)[,1:5]

assayNames(spe)
colData(spe)[1:5,1:5]
rowData(spe)[1:5,1:5]

metaData(spe)$NegProbes[,1:5]

plotSampleInfo(spe, column2plot = c("SlideName","grossRegion", "population", "Group"))

head(spe)

colData(spe)$QCFlags
# No QC flags

# GENE QC

spe <- addPerROIQC(spe, rm_genes = TRUE)
# used defualt parameters -> min_count = 5 and sample_fraction = 0.9


plotGeneQC(spe, ordannots = "grossRegion", col = grossRegion, point_size = 2)

spe <- addPerROIQC(spe)
plotGeneQC(spe)
# not sure whats happening here ngl

# ROI QC
colData(spe)$seqFilt <- colData(spe)$SequencingSaturation < 90 
plotROIQC(spe,  x_axis = "AOISurfaceArea", x_lab = "AOI Surface Area", y_axis = "lib_size", y_lab = "Library size", col = Group)
plotROIQC(spe,  x_axis = "AOISurfaceArea", x_lab = "AOI Surface Area", y_axis = "lib_size", y_lab = "Library size", col = SlideName)
plotROIQC(spe,  x_axis = "AOISurfaceArea", x_lab = "AOI Surface Area", y_axis = "lib_size", y_lab = "Library size", col = grossRegion)
plotROIQC(spe,  x_axis = "AOISurfaceArea", x_lab = "AOI Surface Area", y_axis = "lib_size", y_lab = "Library size", col = population)
plotROIQC(spe,  x_axis = "AOISurfaceArea", x_lab = "AOI Surface Area", y_axis = "lib_size", y_lab = "Library size", col = seqFilt)
# few different labels on the ROIQC plots. Clear separation seen when labelled with population
plotROIQC(spe,  x_axis = "AOISurfaceArea", x_threshold = 6500, x_lab = "AOI Surface Area", y_axis = "lib_size", y_lab = "Library size", col = population)
# whole data set had surface area over 6500. not sure what min is


hist(colData(spe)$AlignedReads)

# Relative log expression distribution

plotRLExpr(spe)
#raw

plotRLExpr(spe, ordannots = "population", assay = 2, color = population)
plotRLExpr(spe, ordannots = "SlideName", assay = 2, color = SlideName)
plotRLExpr(spe, ordannots = "grossRegion", assay = 2, color = grossRegion)
plotRLExpr(spe, ordannots = "Group", assay = 2, color = Group)
# some separation between slides, clear separation between population
# little difference when looking at gross region and group

# Dimension reduction

library(scater)
set.seed(100)

spe <- scater::runPCA(spe)

pca_results <- reducedDim(spe, "PCA")
# 'You're computing too large a percentage of total singular values, use a standard svd instead'

drawPCA(spe, precomputed = pca_results, color = population)
drawPCA(spe, precomputed = pca_results, color = SlideName)
drawPCA(spe, precomputed = pca_results, color = grossRegion)
drawPCA(spe, precomputed = pca_results, color = Group)

plotScreePCA(spe, precomputed = pca_results)

plotPairPCA(spe, col = population, precomputed = pca_results, n_dimension = 4)
plotPairPCA(spe, col = SlideName, precomputed = pca_results, n_dimension = 4)
plotPairPCA(spe, col = grossRegion, precomputed = pca_results, n_dimension = 4)
plotPairPCA(spe, col = Group, precomputed = pca_results, n_dimension = 4)

plotPCAbiplot(spe, n_loadings = 10, precomputed = pca_results, col = grossRegion)
plotPCAbiplot(spe, n_loadings = 10, precomputed = pca_results, col = Group)
plotPCAbiplot(spe, n_loadings = 10, precomputed = pca_results, col = SlideName)
plotPCAbiplot(spe, n_loadings = 10, precomputed = pca_results, col = population)

# MDS

standR::plotMDS(spe, assay = 2, color = grossRegion)
# can make more if needed, does similar thing to pca plots

# UMAP
set.seed(100)

spe <- scater::runUMAP(spe, dimred = "PCA")

plotDR(spe, dimred = "UMAP", col = population)
plotDR(spe, dimred = "UMAP", col = grossRegion)
plotDR(spe, dimred = "UMAP", col = SlideName)
plotDR(spe, dimred = "UMAP", col = Group)

# Normalisation

spe_tmm <- geomxNorm(spe, method = "TMM")
plotRLExpr(spe_tmm, assay = 2, color = population) + ggtitle("TMM")

set.seed(100)

spe_tmm <- scater::runPCA(spe_tmm)
# 'You're computing too large a percentage of total singular values, use a standard svd instead'

pca_results_tmm <- reducedDim(spe_tmm, "PCA")

plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = Group)
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = grossRegion)
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = SlideName)
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = population)

# batch correction

spe <- findNCGs(spe, batch_name = "SlideName", top_n = 300)

for(i in seq(5)){
  spe_ruv <- geomxBatchCorrection(spe, factors = c("Group","grossRegion", "population"), 
                                  NCGs = metadata(spe)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
}

for(i in seq(5)){
  spe_ruv <- geomxBatchCorrection(spe, factors = "Group", 
                                  NCGs = metadata(spe)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
}

for(i in seq(5)){
  spe_ruv <- geomxBatchCorrection(spe, factors = "population", 
                                  NCGs = metadata(spe)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = population, title = paste0("k = ", i)))
  
}

# deconvolution

BiocManager::install("SpatialDecon")
library(SpatialDecon)

norm <- assay(spe_tmm)

duplicates <- c("D830030K20Rik", "Gm10406", "LOC118568634")
subCountData <- countData[-which(countData$TargetName %in% duplicates),]  
rownames(subCountData) <- subCountData$ProbeDisplayName

norm <- rbind(norm,subCountData[grep("NegProbe",subCountData$TargetName),14:108])
negSet<- rownames(subCountData[grep("NegProbe",subCountData$TargetName),14:108])
# row.names(norm)[(nrow(norm) - 209):nrow(norm)] <- paste("NegProbe-WX",1:210)
bg2 = derive_GeoMx_background(norm = norm,
                              probepool = rep(1, nrow(norm)),negnames = negSet)
norm <- as.matrix(assay(spe_tmm))

mousebrain <- download_profile_matrix(species = "Mouse",
                                      age_group = "Adult", 
                                      matrixname = "Brain_AllenBrainAtlas")

heatmap(sweep(mousebrain, 1, apply(mousebrain, 1, max), "/"),
        labRow = NA, margins = c(10, 5))

res <- spatialdecon(norm = norm,
                    bg = bg2,
                    X = mousebrain,
                    align_genes = TRUE)

samples_subset <- colnames(spe_tmm)[colData(spe_tmm)$grossRegion %in%  c("CA1", "EC")]

subset_prop <- res$prop_of_all[,samples_subset]

spe_sub <- spe_tmm[,samples_subset]

BiocManager::install("speckle")
library(speckle)
library(ggrepel)

subset_prop %>%
  as.data.frame() %>%
  rownames_to_column("CellTypes") %>%
  gather(samples, prop, -CellTypes) %>%
  ggplot(aes(samples, prop, fill = CellTypes)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = .7) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")
‹

# differential proportion analysis


BiocManager::install("speckle")
library(speckle)


propslist <- convertDataToList(subset_prop, 
                               data.type = c("proportions"),
                               transform="asin",
                               scale.fac=colData(spe_sub)$AOINucleiCount)

design <- model.matrix(~ 0 + grossRegion + SlideName, data = as.data.frame(colData(spe_sub)))

colnames(design) <- str_remove(colnames(design), pattern = "grossRegion") %>%
  str_replace_all(., " ", "_")

contr <- makeContrasts(CA1 - EC,levels=design)


outs <- propeller.ttest(propslist, design, contr, robust=TRUE,trend=FALSE, sort=TRUE)

diff_ct <- outs %>% 
  filter(FDR < 0.05) %>%
  rownames()

colData(spe_sub)$samples_id <- rownames(colData(spe_sub))

subset_prop[diff_ct,] %>%
  as.data.frame() %>%
  rownames_to_column("CellTypes") %>%
  gather(samples, prop, -CellTypes) %>%
  left_join(as.data.frame(colData(spe_sub)), by = c("samples"="samples_id")) %>%
  ggplot(aes(x = grossRegion, y = prop, fill = population)) +
  geom_violin() +
  facet_wrap(~CellTypes) +
  theme_bw() +
  xlab("gross region") +
  ylab("Proportion")

# comment