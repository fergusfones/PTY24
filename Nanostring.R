# Nanostring script
# Authors - Fergus Fones (University of Exeter), Joshua Harvey (University of Exeter Medical School)



setwd("/Users/fergusfones/Desktop/Nanostring/")
getwd()

#packages
library(standR)
library(SpatialExperiment)
library(ggplot2)
library(ggalluvial)
library(tidyverse)
library(readr)
library(ggpubr)
library(edgeR)
library(limma)
library(GeomxTools)
library(SpatialDecon)
library(fastDummies)  
library(reshape2)
library(dplyr)
library(ggrepel)
library(scater)
library(DT)


metaData <- read.csv("preProcessedMetaData_filtered.csv",header = T)
countDataa <- read.csv("preProcessedCountData_filtered.csv",header = T,)
pathology <- read.csv("PathologyPG5.csv", header = T,)
rownames(pathology) <- pathology$Histology.no.
rownames(metaData) <- metaData$SegmentDisplayName
colnames(countDataa)[14:108] <- gsub("X","",colnames(countDataa)[14:108])

duplicates <- c("D830030K20Rik", "Gm10406", "LOC118568634")
subCountData <- countDataa[-which(countDataa$TargetName %in% duplicates),]  
rownames(subCountData) <- subCountData$ProbeDisplayName

spe <- readGeoMx(countFile = subCountData[,c(4,14:108)],sampleAnnoFile = metaData,
                 featureAnnoFile = subCountData[,-c(14:108)],rmNegProbe = T) 

# general:

plotSampleInfo(spe, column2plot = c("Group", "SlideName", "grossRegion", "population"))

head(spe)

colData(spe)$QCFlags
# No QC flags

# GENE QC

spe <- addPerROIQC(spe, rm_genes = TRUE)
# used defualt parameters -> min_count = 5 and sample_fraction = 0.9
# (filtering genes with low express and below threshold in more than 90% of samples)


plotGeneQC(spe, ordannots = "grossRegion", col = grossRegion, point_size = 2, top_n = 10)
# only 5 genes removed

spe <- addPerROIQC(spe)

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

plotROIQC(spe,  x_axis = "AOINucleiCount",  x_lab = "AOINucleiCount", y_axis = "lib_size", y_lab = "Library size", col = population)

#qc <- colData(spe)$AOINucleiCount > 150
# could do for lib size

#table(qc)

hist(colData(spe)$AlignedReads,
     xlab = "Read Count",
     main = "Aligned Reads",
     xlim = c(0, 3e7), 
     breaks = "FD")

mean((colData(spe)$AlignedReads))

hist(colData(spe)$AOINucleiCount,
     xlab = "Nuclei Count",
     main = "AOI Nuclei Count",
     breaks = "FD",
     xlim = c(0,2000),
     xaxt = 'n')
axis(side=1, at=seq(0, 2000, by=200))

mean((colData(spe)$AOINucleiCount))
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

set.seed(100)

spe <- scater::runPCA(spe)

pca_results <- reducedDim(spe, "PCA")
# 'You're computing too large a percentage of total singular values, use a standard svd instead'

drawPCA(spe, precomputed = pca_results, color = population)
drawPCA(spe, precomputed = pca_results, color = SlideName)
drawPCA(spe, precomputed = pca_results, color = grossRegion)
drawPCA(spe, precomputed = pca_results, color = Group)
drawPCA(spe, precomputed = pca_results, color = Group, dims = c(3,4))

plotScreePCA(spe, precomputed = pca_results)

plotPairPCA(spe, col = population, precomputed = pca_results, n_dimension = 4)
plotPairPCA(spe, col = SlideName, precomputed = pca_results, n_dimension = 4)
plotPairPCA(spe, col = grossRegion, precomputed = pca_results, n_dimension = 4)
plotPairPCA(spe, col = Group, precomputed = pca_results, n_dimension = 4)

plotPCAbiplot(spe, n_loadings = 10, precomputed = pca_results, col = grossRegion)
plotPCAbiplot(spe, n_loadings = 10, precomputed = pca_results, col = Group)
plotPCAbiplot(spe, n_loadings = 10, precomputed = pca_results, col = SlideName)
plotPCAbiplot(spe, n_loadings = 10, precomputed = pca_results, col = population)


#
#
#
#
#
# Normalisation

spe_tmm <- geomxNorm(spe, method = "TMM")
spe_tmm <- addPerROIQC(spe_tmm, rm_genes = FALSE)

plotRLExpr(spe_tmm, assay = 2, color = Group) + ggtitle("TMM")
plotRLExpr(spe_tmm, assay = 2, color = SlideName) + ggtitle("TMM")
plotRLExpr(spe_tmm, assay = 2, color = population) + ggtitle("TMM")
plotRLExpr(spe_tmm, assay = 2, color = grossRegion) + ggtitle("TMM")

set.seed(100)

spe_tmm <- scater::runPCA(spe_tmm)

pca_results_tmm <- reducedDim(spe_tmm, "PCA")


plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = Group)

plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = grossRegion)
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = SlideName)
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = population)

plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = lib_size, n_dimension = 6)


# deconvolution

NormCountData <- assay(spe_tmm)
# NormCountData<- NormCountData[-which(NormCountData$TargetName %in% duplicates),]
# rownames(NormCountData) <- subCountData$ProbeDisplayName
norm <- NormCountData
negSet<- subCountData[grep("NegProbe",subCountData$TargetName),14:108]
negNames <- subCountData[grep("NegProbe",subCountData$TargetName),"ProbeDisplayName"]
rownames(negSet) <- negNames
norm <- rbind(norm,negSet)

bg2 = derive_GeoMx_background(norm = norm,
                              probepool = rep(1, nrow(norm)),negnames = negNames)

mousebrain <- download_profile_matrix(species = "Mouse",
                                      age_group = "Adult", 
                                      matrixname = "Brain_AllenBrainAtlas")
dim(norm)
norm <- norm[-which(rownames(norm) %in% negNames),]
subCountData <- subCountData[which(subCountData$TargetName %in% rownames(norm)),]
rownames(subCountData) <- subCountData$TargetName
raw <- subCountData[,14:108]

resNormRaw <- spatialdecon(norm = as.matrix(norm),
                           raw = as.matrix(raw),
                           bg = bg2,
                           X = mousebrain,
                           align_genes = TRUE)

save(resNormRaw, file = "/Users/fergusfones/Desktop/Nanostring/resNormRaw.Rdata")
load(file = "/Users/fergusfones/Desktop/Nanostring/resNormRaw.Rdata")

cellConvert <- read.csv("CellLabels_Allen.csv",header = T,skip = 1)
cellConvert <- cellConvert[,-1]
cellConvert$Main[which(cellConvert$Main %in% c("Glutamatergic","GABAergic"))] <- "Neuronal"
cellConvert$cellLabs <- gsub("[/-]", ".", cellConvert$Granular)
cellConvert$cellLabs <- gsub(" ", ".", cellConvert$cellLabs)
table(cellConvert$cellLabs %in% rownames(resNormRaw$prop_of_all))
rownames(cellConvert) <- cellConvert$cellLabs
cellSubProp <- as.data.frame(resNormRaw$prop_of_all)
cellSubProp <- cbind(cellSubProp,as.character(cellConvert[rownames(cellSubProp),"Main"]))
colnames(cellSubProp)[ncol(cellSubProp)] <- "BroadCell"
cellSubPropLong<- melt(cellSubProp,id.vars = "BroadCell")
cellSubPropLong <- cellSubPropLong %>% group_by(variable,BroadCell) %>% reframe(sum(value))
cellSubProp  <- dcast(cellSubPropLong,variable~BroadCell)

colData(spe)
colData(spe)[42:45] <- cellSubProp[, 2:5]
colnames(colData(spe))[42:45] <- c( "Astro" , "Immune_Vascular" , "Neuronal" , "Oligo")

#####################################################

# batch correction

# subsetting before RUV corrections

spe_CA1_neun <- spe[,colData(spe)$grossRegion == "CA1" &
                      colData(spe)$population == "neun"]

spe_CA1_rest <- spe[,colData(spe)$grossRegion == "CA1" &
                      colData(spe)$population == "rest"]

spe_EC_neun <- spe[,colData(spe)$grossRegion == "EC" &
                     colData(spe)$population == "neun"]

spe_EC_rest <- spe[,colData(spe)$grossRegion == "EC" &
                     colData(spe)$population == "rest"]

spe_CA1_neun <- findNCGs(spe_CA1_neun, batch_name = "SlideName", top_n = 300)
spe_CA1_rest <- findNCGs(spe_CA1_rest, batch_name = "SlideName", top_n = 300)
spe_EC_neun <- findNCGs(spe_EC_neun, batch_name = "SlideName", top_n = 300)
spe_EC_rest <- findNCGs(spe_EC_rest, batch_name = "SlideName", top_n = 300)

spe <- scater::runPCA(spe)
pca_results <- reducedDim(spe, "PCA")
plotPairPCA(spe, precomputed = pca_results, color = Astro)

for(i in seq(5)){
  spe_CA1_neun_ruv_post <- geomxBatchCorrection(spe_CA1_neun, factors = c("Group", "Astro", "Oligo", "Immune_Vascular"), 
                                                NCGs = metadata(spe_CA1_neun)$NCGs, k = i)
  
  print(plotPairPCA(spe_CA1_neun_ruv_post, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
}

for(i in seq(5)){
  spe_CA1_rest_ruv_post <- geomxBatchCorrection(spe_CA1_rest, factors = c("Group", "Astro", "Oligo", "Immune_Vascular"), 
                                                NCGs = metadata(spe_CA1_rest)$NCGs, k = i)
  
  print(plotPairPCA(spe_CA1_rest_ruv_post, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
}

for(i in seq(5)){
  spe_EC_neun_ruv_post <- geomxBatchCorrection(spe_EC_neun, factors = c("Group", "Astro", "Oligo", "Immune_Vascular"), 
                                               NCGs = metadata(spe_EC_neun)$NCGs, k = i)
  print(plotPairPCA(spe_CA1_neun_ruv_post, assay = 2, n_dimension = 4, color = Neuronal, title = paste0("k = ", i))+scale_fill_viridis_c())
  
  
}

for(i in seq(5)){
  spe_EC_rest_ruv_post <- geomxBatchCorrection(spe_EC_rest, factors = c("Group", "Astro", "Oligo", "Immune_Vascular"), 
                                               NCGs = metadata(spe_EC_rest)$NCGs, k = i)
  
  print(plotPairPCA(spe_EC_rest_ruv_post, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
}


plotPairPCA(spe_CA1_neun_ruv_post, assay = 2, color = Group, title = "spe_CA1_neun_ruv_post", n_dimension = 4)
plotPairPCA(spe_CA1_rest_ruv_post, assay = 2, color = Group, title = "spe_CA1_rest_ruv_post", n_dimension = 4)
plotPairPCA(spe_EC_neun_ruv_post, assay = 2, color = Group, title = "spe_EC_neun_ruv_post", n_dimension = 4)
plotPairPCA(spe_EC_rest_ruv_post, assay = 2, color = Group, title = "spe_EC_rest_ruv_post", n_dimension = 4)



#####################################################
# cor PC test


load("corPCtest.R")
corPCtest(speOb = spe_tmm, pcaOb =  pca_results_tmm, nPC = 10, colLabs = c("SlideName", "lib_size", "AOISurfaceArea", "population", "grossRegion", "Group"))

print(plotPairPCA(spe_tmm, assay = 2, n_dimension = 4, color = lib_size))

corPCtest(speOb = spe_tmm, pcaOb =  pca_results_tmm, nPC = 10, colLabs = c("SlideName", "lib_size", "AOISurfaceArea", "population", "grossRegion", "Group"))


assay(spe_tmm,2)
pca_results_ruv <- reducedDim(assay(spe_ruv,2), "PCA")
corPCtest(speOb = spe_ruv, pcaOb =  pca_results_ruv, nPC = 10, colLabs = c("SlideName", "lib_size", "AOISurfaceArea", "population", "grossRegion", "Group"))

corPCtest(speOb = spe_lrb, pcaOb =  PCA_results_spel, nPC = 10, colLabs = c("SlideName", "lib_size", "AOISurfaceArea", "population", "grossRegion", "Group"))




################################################################





#### Differential gene expression


dge <- SE2DGEList(spe_CA1_neun_ruv_post)


design <- model.matrix(~0 + Group + Astro + Oligo + Immune_Vascular + ruv_W1 + ruv_W2 + ruv_W3 , data = colData(spe_CA1_neun_ruv_post))

#table(colData(spe_ruv)$population,colData(spe_ruv)$grossRegion)

#colnames(design)

contr.matrix <- makeContrasts(
  BvT = GroupCC - GroupWW,
  levels = colnames(design))

#keep <- filterByExpr(dge, design)
#table(keep)

dge <- estimateDisp(dge, design = design, robust = TRUE)

plotBCV(dge, ylim = c(0, 1.3))
bcv_df <- data.frame(
  'BCV' = sqrt(dge$tagwise.dispersion),
  'AveLogCPM' = dge$AveLogCPM,
  'gene_id' = rownames(dge)
)


highbcv <- bcv_df$BCV > 0.6
highbcv_df <- bcv_df[highbcv, ]
points(highbcv_df$AveLogCPM, highbcv_df$BCV, col = "red")
text(highbcv_df$AveLogCPM, highbcv_df$BCV, labels = highbcv_df$gene_id, pos = 4,
     main = "pop bcv")

v <- voom(dge, design)

corfit <- duplicateCorrelation(v, design, block = colData(spe_CA1_rest_ruv_post)$Histology.no.)

v <- voom(dge, design,block = colData(spe_CA1_rest_ruv_post)$Histology.no., correlation =
            corfit$consensus)

corfit <- duplicateCorrelation(v, design, block = colData(spe_CA1_rest_ruv_post)$Histology.no.)

fit <- lmFit(v, design, block = colData(spe_CA1_rest_ruv_post)$Histology.no., correlation =
               corfit$consensus)

fit2 <- contrasts.fit(fit, contr.matrix)

fit2 <- eBayes(fit2)


saveRDS(fit2, file = "/Users/fergusfones/Desktop/Nanostring/fit2_CA1_rest.Rdata")


fit2_CA1_neun <- readRDS(file = "/Users/fergusfones/Desktop/Nanostring/fit2_CA1_neun.Rdata")
fit2_EC_neun <- readRDS(file = "/Users/fergusfones/Desktop/Nanostring/fit2_EC_neun.Rdata")
fit2_CA1_rest <- readRDS(file = "/Users/fergusfones/Desktop/Nanostring/fit2_CA1_rest.Rdata")
fit2_EC_rest <- readRDS(file = "/Users/fergusfones/Desktop/Nanostring/fit2_EC_rest.Rdata")

results_fit2_CA1neun<- decideTests(fit2_CA1_neun)
summary_fit2_CA1neun <- summary(results_fit2_CA1neun)

summary_fit2_CA1neun

results_fit2_ECrest <- decideTests(fit2_CA1_neun)
summary_fit2_ECrest <- summary(results_fit2_ECrest)

summary_fit2_ECrest


de_genes_toptable_CA1 <- topTable(fit2_CA1_neun, coef = 1, sort.by = "p", n = Inf,p.value = 0.05, adjust.method = "fdr",lfc = 0.5) 
de_genes_toptable_EC <- topTable(fit2_EC_neun, coef = 1, sort.by = "p", n = Inf,p.value = 0.05, adjust.method = "fdr",lfc = 0.5) 
de_genes_toptable_CA1_rest <- topTable(fit2_CA1_rest, coef = 1, sort.by = "p", n = Inf,p.value = 0.05,adjust.method = "fdr",lfc = 0.5) 
de_genes_toptable_EC_rest <- topTable(fit2_EC_rest, coef = 1, sort.by = "p", n = Inf,p.value = 0.05, adjust.method = "fdr",lfc = 0.5) 

de_results_CA1 <- topTable(fit2_CA1_neun, coef = 1, sort.by = "P", n = Inf)
de_results_EC <- topTable(fit2_EC_neun, coef = 1, sort.by = "P", n = Inf)
de_results_CA1_rest <- topTable(fit2_CA1_rest, coef = 1, sort.by = "P", n = Inf)
de_results_EC_rest <- topTable(fit2_EC_rest, coef = 1, sort.by = "P", n = Inf)


geneIndex <- rownames(de_results_CA1[which(de_results_CA1$adj.P.Val < 0.05 & de_results_CA1$P.Value < 0.05 &  abs(de_results_CA1$logFC) > 0.5),])
geneIndex2 <- rownames(de_results_EC_rest[which(de_results_EC_rest$adj.P.Val < 0.05 & de_results_EC_rest$P.Value < 0.05 &  abs(de_results_EC_rest$logFC) > 0.5),])
geneIndexx <- rownames(de_results_CA1)

plot(de_results_CA1[geneIndex,"logFC"],de_results_EC[geneIndex,"logFC"]) 
text(de_results_CA1[geneIndex, "logFC"], de_results_EC[geneIndex, "logFC"], 
     labels = geneIndex,
     pos = 4, cex = 0.8, col = "red")
abline(h= 0)
abline(v= 0)

cor.test(de_results_CA1[geneIndex,"logFC"],de_results_CA1_rest[geneIndex,"logFC"])
binom.test(146+80,nrow(de_results_CA1[geneIndex,]))
table(sign(de_results_CA1[geneIndex,"logFC"]),sign(de_results_CA1_rest[geneIndex,"logFC"]))

cor.test(de_results_CA1[geneIndex,"logFC"],de_results_EC[geneIndex,"logFC"])
binom.test(114+78,nrow(de_results_CA1[geneIndex,]))
table(sign(de_results_CA1[geneIndex,"logFC"]),sign(de_results_EC[geneIndex,"logFC"]))

cor.test(de_results_EC[geneIndex,"logFC"],de_results_EC_rest[geneIndex,"logFC"])
binom.test(99+85,nrow(de_results_EC[geneIndex,]))
table(sign(de_results_EC[geneIndex,"logFC"]),sign(de_results_EC_rest[geneIndex,"logFC"]))

cor.test(de_results_CA1_rest[geneIndex,"logFC"],de_results_EC_rest[geneIndex,"logFC"])
binom.test(112+68,nrow(de_results_CA1_rest[geneIndex,]))
table(sign(de_results_CA1_rest[geneIndex,"logFC"]),sign(de_results_EC_rest[geneIndex,"logFC"]))


gh <- data.frame(de_results_CA1[geneIndex,"logFC"],
                 de_results_EC[geneIndex,"logFC"],
                 row.names = rownames(geneIndex))

hj <- data.frame(geneIndex, rowSums(abs(gh)))
ll <- data.frame(geneIndex, rowSums(gh))


hj <- dplyr::arrange(hj, desc(rowSums.abs.gh..))
ll <- dplyr::arrange(ll, desc(rowSums.gh.))


hj[1:20,]
ll[c(1:10, 283:292),]


de_results_EC_rest %>% 
  mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                     ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
  ggplot(aes(logFC, -log10(P.Value), col = DE)) + 
  geom_point(shape = 1, size = 1) + 
  geom_text_repel(data = de_genes_toptable_EC_rest %>% 
                    mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                                       ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
                    rownames_to_column(), aes(label = rowname)) +
  theme_bw() +
  xlab("Log-fold-change") +
  ylab("-log10 P value") +
  ggtitle("WW v CC (spe_EC_rest_ruv_post)") +
  scale_color_manual(values = c("blue","gray","red")) +
  theme(text = element_text(size=15))


updn_cols <- c(RColorBrewer::brewer.pal(6, 'Greens')[2], RColorBrewer::brewer.pal(6, 'Purples')[2])

de_genes_toptable_CA1 %>% 
  dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
  DT::datatable(caption = ' WW v CC de_genes_toptable_CA1') %>%
  DT::formatStyle('logFC',
                  valueColumns = 'logFC',
                  backgroundColor = DT::styleInterval(0, rev(updn_cols))) %>%
  DT::formatSignif(1:4, digits = 4)

de_genes_toptable_EC %>% 
  dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
  DT::datatable(caption = ' WW v CC de_genes_toptable_EC') %>%
  DT::formatStyle('logFC',
                  valueColumns = 'logFC',
                  backgroundColor = DT::styleInterval(0, rev(updn_cols))) %>%
  DT::formatSignif(1:4, digits = 4)





######################
# Target gene screening



de_results_CA1 <- topTable(fit2_CA1_neun, coef = 1, sort.by = "P", n = Inf)
de_results_EC <- topTable(fit2_EC_neun, coef = 1, sort.by = "P", n = Inf)

colnames(de_results_CA1) <- paste(colnames(de_results_CA1), "CA1",sep = "_")
CA1_EC <- cbind(de_results_CA1, de_results_EC[rownames(de_results_CA1), ])
sign(CA1_EC$logFC)  == sign(CA1_EC$logFC_CA1)

CA1_EC$sign <- sign(CA1_EC$logFC)  == sign(CA1_EC$logFC_CA1)

sig_FDR <- CA1_EC[which(CA1_EC$adj.P.Val_CA1 < 0.05 & abs(CA1_EC$logFC_CA1) > 0.5),]

artifacts <- c("Wdr60", "Esyt2", "Ncapg2", "Ptprn2", "Fgf14")
ff<- CA1_EC[which(CA1_EC$adj.P.Val_CA1 < 0.05 & abs(CA1_EC$logFC_CA1) > 0.5 & CA1_EC$adj.P.Val < 0.05 & abs(CA1_EC$logFC) > 0.5),]

sig_FDR <- sig_FDR[which(sig_FDR$sign == TRUE),]

ff %>% 
  mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                     ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
  ggplot(aes(logFC, -log10(P.Value), col = DE)) + 
  geom_point(shape = 1, size = 1) + 
  geom_text_repel(data = ff %>% 
                    mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                                       ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
                    rownames_to_column(), aes(label = rowname)) +
  theme_bw() +
  xlab("Log-fold-change") +
  ylab("-log10 P value") +
  ggtitle("significant genes") +
  scale_color_manual(values = c("blue","gray","red")) +
  theme(text = element_text(size=15))

sigEC <- de_results_EC[which(de_results_EC$adj.P.Val < 0.05 & abs(de_results_EC$logFC) > 0.5 ),]

de_results_EC[de_results_CA1$logFC < -0.5 & de_results_EC$logFC > 0.5 & de_results_EC$adj.P.Val <0.05 ,]

