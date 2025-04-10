# The code performs a WGCNA analysis
# 
# Made using this tutorial : https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0
# And with help from Giulia Pegoraro, University of Exeter Medical School and Lachlan McBean, University of Exeter Medical School

###########################################
setwd("/Users/fergusfones/Desktop/Nanostring/")
getwd()


#install.packages("BiocManager")
#BiocManager::install(c("WGCNA", "devtools"))
library(devtools)
install_github('jdrudolph/PerseusR')
library(tidyverse)
library(WGCNA)
library(PerseusR)
library(DESeq2)
library(dplyr)
library(pheatmap)
#install.packages("RColorBrewer")
library(RColorBrewer)
library(magrittr)
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(purrr)
#iocManager::install("VennDiagram")
library(VennDiagram)
###########################################



metaData <- read.csv("preProcessedMetaData_filtered.csv",header = T)
countDataa <- read.csv("preProcessedCountData_filtered.csv",header = T,)
pathology <- read.csv("PathologyPG5.csv", header = T,)
rownames(pathology) <- pathology$Histology.no.
rownames(metaData) <- metaData$SegmentDisplayName
colnames(countDataa)[14:108] <- gsub("X","",colnames(countDataa)[14:108])

duplicates <- c("D830030K20Rik", "Gm10406", "LOC118568634")
subCountData <- countDataa[-which(countDataa$TargetName %in% duplicates),]  
rownames(subCountData) <- subCountData$ProbeDisplayName



# need to transpose data, so gene probes are in the columns and 'treatments' are in the rows, in a count matrix
WGCNA_matrix <- subCountData[,c(4,14:108)]

WGCNA_matrix <- WGCNA_matrix %>% 
  filter(., TargetName != "NegProbe-WTX")

rownames(WGCNA_matrix) <- WGCNA_matrix[, 1]

WGCNA_matrix <- WGCNA_matrix[, 2:96]

WGCNA_matrix <- as.matrix(WGCNA_matrix)

dim(WGCNA_matrix)
# 19959 genes, 95 samples



#Subsetting before normalisation


WGCNA_matrix_CA1_neun <- WGCNA_matrix[, grepl("neun", colnames(WGCNA_matrix)) & !grepl("rest", colnames(WGCNA_matrix)) & grepl("CA1", colnames(WGCNA_matrix)) | colnames(WGCNA_matrix) == "TargetName"]


WGCNA_matrix_CA1_rest <- WGCNA_matrix[, grepl("rest", colnames(WGCNA_matrix)) & !grepl("neun", colnames(WGCNA_matrix)) & grepl("CA1", colnames(WGCNA_matrix)) | colnames(WGCNA_matrix) == "TargetName"]


WGCNA_matrix_EC_neun <- WGCNA_matrix[, grepl("neun", colnames(WGCNA_matrix)) & !grepl("rest", colnames(WGCNA_matrix)) & grepl("EC", colnames(WGCNA_matrix)) | colnames(WGCNA_matrix) == "TargetName"]


WGCNA_matrix_EC_rest <- WGCNA_matrix[, grepl("rest", colnames(WGCNA_matrix)) & !grepl("neun", colnames(WGCNA_matrix)) & grepl("EC", colnames(WGCNA_matrix)) | colnames(WGCNA_matrix) == "TargetName"]

#
#
#
#
#

#                  Normalising data with DESeq
###

##
#
#
#
#
#
#

meta_df <- read.csv("preProcessedMetaData_filtered.csv",header = T)
meta_df <- meta_df[, c(6, 32, 35, 36, 38 )]


rownames(meta_df) <- meta_df[,1]

# sorting for CA1 neun
meta_df_CA1_neun <- meta_df %>%
  filter(., population == "neun" & grepl("CA1", rownames(meta_df)))
  
# sorting for CA1 rest
meta_df_CA1_rest <- meta_df %>%
  filter(., population == "rest" & grepl("CA1", rownames(meta_df)))

# sorting for EC neun
meta_df_EC_neun <- meta_df %>%
  filter(., population == "neun" & grepl("EC", rownames(meta_df)))

# sorting for EC rest
meta_df_EC_rest <- meta_df %>%
  filter(., population == "rest" & grepl("EC", rownames(meta_df)))


dds <- DESeqDataSetFromMatrix(round(WGCNA_matrix_CA1_neun),
                              meta_df_CA1_neun,
                              design = ~Group)
dds <- DESeq(dds)

# DEseq pipeline ends


vsd <- varianceStabilizingTransformation(dds)
# retains metadata, can be used for PCA, visualisation, further processing
# linear regression fit was substituted for a local regression fit by the function

wpn_vsd <- getVarianceStabilizedData(dds)

rv_wpn <- rowVars(wpn_vsd) 
summary(rv_wpn)


q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= 95 quantile reduces dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]

dim(expr_normalized)
#998, 95
# for Ca1 nuen 998,24
# for ec neun 998,23
# for Ca1 rest 998,24
# for ec rest 998,24


expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 95 quantile Expression_EC_rest",
    x = "treatment",
    y = "normalized expression"
  )


  
  
  #
  #
  ##
  #
  #
  #
  #
  #
  # correlation matrix
  
  
  cor_matrix <-  t(expr_normalized)
  cor_matrix <- cor(cor_matrix, method = "spearman")
  
  colours <- colorRampPalette(brewer.pal(9, "Blues"))(225)
  
  pheatmap(cor_matrix, col = colours, main = "gene expression correlation", show_colnames = F, show_rownames = F)
  
  
  
  # hierarchical clustering 
  
htree <- hclust(dist(t(expr_normalized)))
plot(htree, xlab = "samples", main = "hierarchical clustering ")


  #
  #
  ##
  #
  #
  #
  #
  #
  ##
  #
  #
  #
  #
  #
  
  #                    starting WGCNA
  ###
  
  ##
  #
  #
  #
  #
  #
  #
# 

input_mat <- t(expr_normalized)

# Allowing multi-threading
allowWGCNAThreads() 

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)


par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red",
)

abline(h = 0.8, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

# try using the estimate power function

picked_power <- 9

temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

cor <- temp_cor




# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

# netwk$colors[netwk$blockGenes[[1]]]
# table(netwk$colors)






##
#
#
#
#
#
#
#          Relate Modules assignments to treatment groups
##
#
#
#
#
#
#
# 


module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]


#write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")


# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)


# Adding other variables
MEs0 <- cbind(MEs0, meta_df_CA1_neun)

### Remove grey module (junk)
MEs0 <- MEs0[,-which(colnames(MEs0)=="MEgrey")]

traits <- meta_df_CA1_neun %>%
  select(- c("SegmentDisplayName", "Histology.no.", "population", "grossRegion"))

#encoding traits
traits_numeric <- ifelse(traits == "WW", 0 ,1) 

# module trait correlation
# adjust rows in MEs0 to just show modules
moduleTraitCor = stats::cor(MEs0[,1:11],traits_numeric , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(MEs0))


textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

# Adjust margins and font sizes
par(mfrow = c(1, 1))
par(mar = c(2, 2, 3, 3), cex.main = 1.5, cex.axis = 1, cex.lab = 0.7)  # Decrease margins and adjust font sizes

# Plot the heatmap

#pheatmap(moduleTraitCor, col = colours, main = "module trait correlation")


labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traits_numeric),
  yLabels = names(MEs0[,1:11]),
  ySymbols = names(MEs0[,1:11]),
  colorLabels = FALSE,
  colors = colorRampPalette(c("blue", "white", "red"))(50),  # Change color scheme if needed
  textMatrix = textMatrix,
  setStdMargins = TRUE,
  plotLegend = TRUE,
  cex.text = 0.6,  # Reduce label size
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)



##
#
#
#
#
#
#
#          Examine expression profiles

##
#
#
#
#
#
#
# 


# how many genes in each module
table(module_df$colors)

# pick out a few modules of interest here
modules_of_interest = c("turquoise")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
expr_normalized[1:5,1:10]

subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")



##
#
#
#
#
#
#
#          Generate and export networks

##
#
#
#
#
#
#
# 





genes_of_interest = module_df %>%
  subset(colors == "turquoise")


expr_of_interest = expr_normalized[genes_of_interest$gene_id,]
#expr_of_interest[1:5,1:5]

TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)

# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "/Users/fergusfones/Desktop/Nanostring/turq_list.tsv",
            delim = "\t")

#
#
#
#
#
##
##
#
#
#
##
#
##
#
# Gene ontology enrichment
##
##
#
#
#
##
#
##
#


# locate which module contains gene of interest and investigate ontology of module
module_df[module_df$gene_id == "Gfap",]
module_df[module_df$gene_id == "Pfkm",]


WGCNA_turq <-module_df$gene_id[module_df$colors == "turquoise"]
WGCNA_brown <-module_df$gene_id[module_df$colors == "brown"]
WGCNA_black <-module_df$gene_id[module_df$colors == "black"]
WGCNA_green <-module_df$gene_id[module_df$colors == "green"]
WGCNA_blue <-module_df$gene_id[module_df$colors == "blue"]
WGCNA_yellow <-module_df$gene_id[module_df$colors == "yellow"]
WGCNA_greenyellow <-module_df$gene_id[module_df$colors == "greenyellow"]
WGCNA_magenta <-module_df$gene_id[module_df$colors == "magenta"]
WGCNA_red <-module_df$gene_id[module_df$colors == "red"]
WGCNA_pink <-module_df$gene_id[module_df$colors == "pink"]
WGCNA_purple <-module_df$gene_id[module_df$colors == "purple"]

library(enrichR)

dbs <- listEnrichrDbs()

head(dbs)


dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023")

module_enriched <- enrichr(WGCNA_yellow, dbs)

# maybe add p value threshold at 0.05

plotEnrich(module_enriched[[2]], showTerms = 20, numChar = 120, y = "Count", orderBy = "P.value")







##
##
#
#
#
##
#
##
#
#          EWCE of module
##
##
#
#
#
##
#
##
#





library(EWCE)
library(ewceData)

mouseCTD <- ewceData::ctd()


#finally subset the list of names for your DEGs
topGenes <- WGCNA_turq


# This tests for enrichment in the mouse single cell dataset
# Start out with 100 for quick running but when we run the analysis "finally" make this 10,000 reps
results <- EWCE::bootstrap_enrichment_test(sct_data = mouseCTD,
                                           sctSpecies = "mouse",
                                           genelistSpecies = "mouse",
                                           hits = topGenes, 
                                           reps = 100,
                                           annotLevel = 1)


#Take out the enrichment results
ewceRes <- results$results


# This method only tests for positive enrichment (sd_from_mean < 0 == NA)
ewceRes[which(ewceRes$sd_from_mean < 0),"sd_from_mean"] <- NA

#annotate significance
ewceRes$sigAnnot <- NA
ewceRes$sigAnnot[ewceRes$q < 0.05] <- "*"
ewceRes$sigAnnot[ewceRes$q < 0.0001] <- "**"
ewceRes$sigAnnot[ewceRes$q < 1e-10] <- "***"

#Finally, plot your cell types
ewceRes %>%
  ggplot(aes(x = as.factor(CellType), y = sd_from_mean))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(ewceRes$sd_from_mean,na.rm = T) + max(ewceRes$sd_from_mean,na.rm = T)/10))+
  geom_text(aes(label = sigAnnot, y = sd_from_mean + max(ewceRes$sd_from_mean,na.rm = T)/20), size = 8)+
  ylab("SD From Mean")+
  xlab("Cell type")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

