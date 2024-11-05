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
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("GeomxTools")
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
# (filtering genes with low express and below threshold in more than 90% of samples)


plotGeneQC(spe, ordannots = "grossRegion", col = grossRegion, point_size = 2, top_n = 10)
# only 5 genes removed

spe <- addPerROIQC(spe)
plotGeneQC(spe)
# i think its adding the removed gene function to spe.

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

standR::plotMDS(spe, assay = 2, color = population)
#can make more if needed, does similar thing to pca plots

# UMAP
set.seed(100)

spe <- scater::runUMAP(spe, dimred = "PCA")

plotDR(spe, dimred = "UMAP", col = population)
plotDR(spe, dimred = "UMAP", col = Group)


# Normalisation

spe_tmm <- geomxNorm(spe, method = "TMM")
plotRLExpr(spe_tmm, assay = 2, color = Group) + ggtitle("TMM")

#spe_cpm <- geomxNorm(spe, method = "CPM")
#plotRLExpr(spe_cpm, assay = 2, color = population) + ggtitle("CPM")

set.seed(100)

spe_tmm <- scater::runPCA(spe_tmm)

#spe_cpm <- scater::runPCA(spe_cpm)
# 'You're computing too large a percentage of total singular values, use a standard svd instead'

pca_results_tmm <- reducedDim(spe_tmm, "PCA")
#pca_results_cpm <- reducedDim(spe_cpm, "PCA")

plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = Group)
plotPairPCA(spe_cpm, precomputed = pca_results_cpm, color = Group)

plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = grossRegion)
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = SlideName)
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = population)

plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = lib_size, n_dimension = 6)

# batch correction

# substetting before RUV corrections

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

for(i in seq(5)){
  spe_CA1_neun_ruv_post <- geomxBatchCorrection(spe_CA1_neun, factors = c("Group","grossRegion"), 
                                                NCGs = metadata(spe_CA1_neun)$NCGs, k = i)
  
  print(plotPairPCA(spe_CA1_neun_ruv_post, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
}

for(i in seq(5)){
  spe_CA1_rest_ruv_post <- geomxBatchCorrection(spe_CA1_rest, factors = c("Group","grossRegion"), 
                                                NCGs = metadata(spe_CA1_rest)$NCGs, k = i)
  
  print(plotPairPCA(spe_CA1_rest_ruv_post, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
}

for(i in seq(5)){
  spe_EC_neun_ruv_post <- geomxBatchCorrection(spe_EC_neun, factors = c("Group","grossRegion"), 
                                               NCGs = metadata(spe_EC_neun)$NCGs, k = i)
  
  print(plotPairPCA(spe_EC_neun_ruv_post, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
}

for(i in seq(5)){
  spe_EC_rest_ruv_post <- geomxBatchCorrection(spe_EC_rest, factors = c("Group","grossRegion"), 
                                               NCGs = metadata(spe_EC_rest)$NCGs, k = i)
  
  print(plotPairPCA(spe_EC_rest_ruv_post, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
}


plotPairPCA(spe_CA1_neun_ruv_post, assay = 2, color = Group, title = "spe_CA1_neun_ruv_post", n_dimension = 4)
plotPairPCA(spe_CA1_rest_ruv_post, assay = 2, color = Group, title = "spe_CA1_rest_ruv_post", n_dimension = 4)
plotPairPCA(spe_EC_neun_ruv_post, assay = 2, color = Group, title = "spe_EC_neun_ruv_post", n_dimension = 4)
plotPairPCA(spe_EC_rest_ruv_post, assay = 2, color = Group, title = "spe_EC_rest_ruv_post", n_dimension = 4)



# deconvolution

BiocManager::install("SpatialDecon")
library(SpatialDecon)

norm <- assay(spe_cpm)

duplicates <- c("D830030K20Rik", "Gm10406", "LOC118568634")
subCountData <- countData[-which(countData$TargetName %in% duplicates),]  
rownames(subCountData) <- subCountData$ProbeDisplayName

norm <- rbind(norm,subCountData[grep("NegProbe",subCountData$TargetName),14:108])
negSet<- rownames(subCountData[grep("NegProbe",subCountData$TargetName),14:108])
# row.names(norm)[(nrow(norm) - 209):nrow(norm)] <- paste("NegProbe-WX",1:210)
bg2 = derive_GeoMx_background(norm = norm,
                              probepool = rep(1, nrow(norm)),negnames = negSet)
norm <- as.matrix(assay(spe_cpm))

mousebrain <- download_profile_matrix(species = "Mouse",
                                      age_group = "Adult", 
                                      matrixname = "Brain_AllenBrainAtlas")

mousebrain2 <- download_profile_matrix(species = "Mouse",
                                       age_group = "Adult", 
                                       matrixname = "Brain_MCA")

#mousebrain3 <- download_profile_matrix(species = "Mouse",
                                      age_group = "Adult", 
                                      matrixname = "Allen_Brain_Atlas_10x_scRNA_2021")


res <- spatialdecon(norm = norm,
                    bg = bg2,
                    X = mousebrain,
                    align_genes = TRUE)

res2 <- spatialdecon(norm = norm,
                    bg = bg2,
                    X = mousebrain2,
                    align_genes = TRUE)

res3 <- spatialdecon(norm = norm,
                     bg = bg2,
                     X = custom_mtx_seurat,
                     align_genes = TRUE)


save(res, file = "/Users/fergusfones/Desktop/Nanostring/res.Rdata")
save(res2, file = "/Users/fergusfones/Desktop/Nanostring/res2.Rdata")
save(res3, file = "/Users/fergusfones/Desktop/Nanostring/res3.Rdata")

load(file = "/Users/fergusfones/Desktop/Nanostring/res.Rdata")
load(file = "/Users/fergusfones/Desktop/Nanostring/res2.Rdata")
load(file = "/Users/fergusfones/Desktop/Nanostring/res3.Rdata")


#colSums(res$prop_of_all)

#samples_subset <- colnames(spe_tmm)[colData(spe_tmm)$grossRegion %in%  c("CA1", "EC")]

#subset_prop <- res$prop_of_all

#spe_sub <- spe_tmm[,samples_subset]

#long_subset_prop <- subset_prop %>%
  #as.data.frame() %>%
  #rownames_to_column("CellTypes") %>%
  #gather(samples, prop, -CellTypes)

BiocManager::install("speckle")
library(speckle)
library(ggrepel)

# restructing data from wide to long
# setup for res 1
res_prop.df <- as.data.frame(res$prop_of_all)

long_res <- res_prop.df %>% 
  rownames_to_column(var = "cell_type") %>%
  pivot_longer(cols = "4_1_CA1_neun":"1_6_EC_rest",
               names_to = "region", 
               values_to = "proportion")

region<- colnames(spe_tmm)

pop_group <- cbind(region, spe_tmm$grossRegion, spe_tmm$population, spe_tmm$Group)

long_ress_added <- left_join(long_res, pop_group, by = 'region', copy = T)
colnames(long_ress_added)[4]<-paste("grossRegion")
colnames(long_ress_added)[5]<-paste("population")
colnames(long_ress_added)[6]<-paste("Group") 

long_ress_added

# setup for res 2

res2_prop.df <- as.data.frame(res2$prop_of_all)

long_res2 <- res2_prop.df %>% 
  rownames_to_column(var = "cell_type") %>%
  pivot_longer(cols = "4_1_CA1_neun":"1_6_EC_rest",
               names_to = "region", 
               values_to = "proportion")

region<- colnames(spe_tmm)

pop_group <- cbind(region, spe_tmm$grossRegion, spe_tmm$population, spe_tmm$Group)

long_ress_added2 <- left_join(long_res2, pop_group, by = 'region', copy = T)
colnames(long_ress_added2)[4]<-paste("grossRegion")
colnames(long_ress_added2)[5]<-paste("population")
colnames(long_ress_added2)[6]<-paste("Group") 

long_ress_added2

#deconvolution of res and res2
 
long_ress_added2 %>%
  ggplot(aes(x = region, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = .7) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")

long_ress_added2 %>%
  filter(grossRegion == "EC") %>%
  ggplot(aes(Group, proportion, fill = population)) +
  geom_violin() +
  facet_wrap(~cell_type) +
  theme_bw() +
  xlab("") +
  ylab("Proportion")

long_ress_added2 %>%
filter(population == "neun") %>%
  ggplot(aes(x = cell_type, y = region, fill = proportion)) +
  geom_tile()+
  theme_bw() +
  xlab("Cell Type") +
  ylab("Individidual/region")+
  scale_fill_viridis_c()

long_ress_added %>%
  ggplot(aes(x = region, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = .7) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")

long_ress_added %>%
  filter(grossRegion == "CA1") %>%
  ggplot(aes(Group, proportion, fill = population)) +
  geom_violin() +
  facet_wrap(~cell_type) +
  theme_bw() +
  xlab("") +
  ylab("Proportion")



# deconvolution using long read ref data


long_res3 <- res3$prop_of_all %>% 
  as.data.frame() %>%
  rownames_to_column(var = "cell_type") %>%
  pivot_longer(cols = "4_1_CA1_neun":"1_6_EC_rest",
               names_to = "region", 
               values_to = "proportion")

long_res3 %>%
  ggplot(aes(x = region, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = .7) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")

ggplot(long_res3[-grep("neun",long_res3$region),],aes(x = region, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = .7) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")

long_res3 %>%
  ggplot(aes(region, proportion, fill = cell_type)) +
  geom_violin() +
  facet_wrap(~cell_type) +
  theme_bw() +
  xlab("") +
  ylab("Proportion")


heatmap(res3$beta, cexCol = 0.5, cexRow = 0.7, margins = c(10,7))
image()


# cor PC test


load("corPCtest.R")
corPCtest(speOb = spe_tmm, pcaOb =  pca_results_tmm, nPC = 10, colLabs = c("SlideName", "lib_size", "AOISurfaceArea", "population", "grossRegion", "Group"))
library(ggplot2)
library(tidyverse)
library(dplyr)
library(fastDummies)  
library(reshape2)

print(plotPairPCA(spe_tmm, assay = 2, n_dimension = 4, color = lib_size))

corPCtest(speOb = spe_tmm, pcaOb =  pca_results_tmm, nPC = 10, colLabs = c("SlideName", "lib_size", "AOISurfaceArea", "population", "grossRegion", "Group"))


assay(spe_tmm,2)
pca_results_ruv <- reducedDim(assay(spe_ruv,2), "PCA")
corPCtest(speOb = spe_ruv, pcaOb =  pca_results_ruv, nPC = 10, colLabs = c("SlideName", "lib_size", "AOISurfaceArea", "population", "grossRegion", "Group"))

corPCtest(speOb = spe_lrb, pcaOb =  PCA_results_spel, nPC = 10, colLabs = c("SlideName", "lib_size", "AOISurfaceArea", "population", "grossRegion", "Group"))




################################################################





#### Differential expression

library(edgeR)
library(limma)


dge <- SE2DGEList(spe_EC_rest_ruv_post)


design <- model.matrix(~0 + Group + ruv_W1 + ruv_W2 + ruv_W3 , data = colData(spe_EC_rest_ruv_post))

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

corfit <- duplicateCorrelation(v, design, block = colData(spe_EC_rest_ruv_post)$Histology.no.)

v <- voom(dge, design,block = colData(spe_EC_rest_ruv_post)$Histology.no., correlation =
            corfit$consensus)

corfit <- duplicateCorrelation(v, design, block = colData(spe_EC_rest_ruv_post)$Histology.no.)

fit <- lmFit(v, design, block = colData(spe_EC_rest_ruv_post)$Histology.no., correlation =
               corfit$consensus)

fit2 <- contrasts.fit(fit, contr.matrix)

fit2 <- eBayes(fit2)

results_fit2<- decideTests(fit2)
summary_fit2 <- summary(results_fit2)

summary_fit2

library(ggrepel)
library(tidyverse)

de_genes_toptable <- topTable(fit2, coef = 1, sort.by = "p", n = Inf,p.value = 0.05,adjust.method = "fdr",lfc = 0.5) 

de_results <- topTable(fit2, coef = 1, sort.by = "P", n = Inf)
# there is no difference between these two objects, but ill keep both.
de_results %>% 
  mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                     ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
  ggplot(aes(AveExpr, logFC, col = DE)) + 
  geom_point(shape = 1, size = 1) + 
  geom_text_repel(data = de_genes_toptable %>% 
                    mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                                       ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
                    rownames_to_column(), aes(label = rowname)) +
  theme_bw() +
  xlab("Average log-expression") +
  ylab("Log-fold-change") +
  ggtitle("WW v CC (spe_CA1_rest_ruv_post)") +
  scale_color_manual(values = c("blue","gray","red")) +
  theme(text = element_text(size=15))

library(DT)

updn_cols <- c(RColorBrewer::brewer.pal(6, 'Greens')[2], RColorBrewer::brewer.pal(6, 'Purples')[2])

de_genes_toptable %>% 
  dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
  DT::datatable(caption = ' WW v CC (limma-voom) spe_EC_neun_ruv_post') %>%
  DT::formatStyle('logFC',
                  valueColumns = 'logFC',
                  backgroundColor = DT::styleInterval(0, rev(updn_cols))) %>%
  DT::formatSignif(1:4, digits = 4)

