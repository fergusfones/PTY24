setwd("/Users/fergusfones/Desktop/Nanostring/")
getwd()
#file_path <- "/Users/fergusfones/Desktop/Nanostring/GSM7021713_C_3m_filtered_feature_bc_matrix.h5"
#h5ls(file_path)
#h5read(file_path, /matrix/data )


#expr_matrix <- h5read(file_path, "/matrix/data")

#genes <- h5read(file_path, "/matrix/features")

#barcodes <- h5read(file_path, "/matrix/barcodes")


#indptr <- h5read(file_path, "/matrix/indptr")

#install.packages("hdf5r")
#library(hdf5r)


read.table("/Users/fergusfones/Desktop/Nanostring/GSE224398_series_matrix.txt")

################### seurat

install.packages('Seurat')
library(Seurat)
library(SpatialDecon)
library(dplyr)


setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))

WTref<- Read10X_h5("/Users/fergusfones/Desktop/Nanostring/GSM7021713_C_3m_filtered_feature_bc_matrix.h5")
TGref<- Read10X_h5("/Users/fergusfones/Desktop/Nanostring/GSM7021716_A_3m_filtered_feature_bc_matrix.h5")

combinedref <- RowMergeSparseMatrices(WTref, TGref)





################## szi kay stuff

long_reads<- read.delim("/Users/fergusfones/Desktop/Nanostring/rTg4510SCN_collapsed_RulesFilter_result_classification_counts.txt")
phenotypes <- read.csv("/Users/fergusfones/Desktop/Nanostring/SCNPhenotype.csv")

Gfap <- long_reads[long_reads$associated_gene == "Gfap",] %>%
  select(contains("mapped")) %>%
  summarise(across(1:16, ~ sum(.x, na.rm = TRUE)))

Gfap %>% reshape2::melt(variable.name = "sample",value.name = "reads") %>% 
  mutate(sample_id = stringr::word(sample,c(1), sep = fixed("_"))) %>% 
  merge(., phenotypes, by = "sample") %>% 
  mutate(group = factor(group, levels = c("WT","TG"))) %>%
  ggplot(., aes(x = group, y = reads, colour = cell)) + geom_point(size = 3) +
  facet_grid(~age)

Gfap %>% reshape2::melt(variable.name = "sample",value.name = "reads") %>% 
  mutate(sample_id = stringr::word(sample,c(1), sep = fixed("_"))) %>% 
  merge(., phenotypes, by = "sample") %>% 
  mutate(group = factor(group, levels = c("WT","TG"))) %>%
  ggplot(., aes(x = group, y = reads, colour = cell)) + geom_boxplot()



library(dplyr)
# n_distinct(long_reads$associated_gene)
# dim(long_reads)
# ncol(long_reads)
# colnames(long_reads)

library(edgeR)
library(ggplot2)

subset_long_reads <- long_reads %>% select(7, 58:68, 70:74)

summed_data <- subset_long_reads %>%
  group_by(associated_gene = subset_long_reads[[1]]) %>%  
  summarise(across(1:16, ~ sum(.x, na.rm = TRUE)))
  

summed_data <- as.data.frame(summed_data)
rownames(summed_data) <- summed_data$associated_gene
summed_data_temp <- summed_data[,-1]

summed_data_cpm <- cpm(summed_data_temp, method = "CPM")

custom_mtx_seurat <- create_profile_matrix(mtx = summed_data_cpm ,
                                           cellAnnots = phenotypes, 
                                           cellTypeCol = "cell", 
                                           cellNameCol = "sample", 
                                           matrixName = "custom_mini_fans",
                                           outDir = NULL, 
                                           normalize = FALSE, 
                                           minCellNum = 5, 
                                           minGenes = 10)


rownames(phenotypes) <- phenotypes$sample
phenotypes[colnames(summed_data),]
