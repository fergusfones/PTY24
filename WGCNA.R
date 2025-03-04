setwd("/Users/fergusfones/Desktop/Nanostring/")
getwd()
# WGCNA

install.packages("BiocManager")
BiocManager::install(c("WGCNA", "devtools"))
library(devtools)
install_github('jdrudolph/PerseusR')
library(tidyverse)
library(WGCNA)
library(PerseusR)
library(DESeq2)
library(dplyr)
library(pheatmap)
install.packages("RColorBrewer")
library(RColorBrewer)

metaData <- read.csv("preProcessedMetaData_filtered.csv",header = T)
countDataa <- read.csv("preProcessedCountData_filtered.csv",header = T,)
pathology <- read.csv("PathologyPG5.csv", header = T,)
rownames(pathology) <- pathology$Histology.no.
rownames(metaData) <- metaData$SegmentDisplayName
colnames(countDataa)[14:108] <- gsub("X","",colnames(countDataa)[14:108])

duplicates <- c("D830030K20Rik", "Gm10406", "LOC118568634")
subCountData <- countDataa[-which(countDataa$TargetName %in% duplicates),]  
rownames(subCountData) <- subCountData$ProbeDisplayName



# do i have to control for covariates like cell types and correct the data like we did in the GeoMX workflow?
# Currently just using CA1 neun samples, current n=24


# need to transpose data, so gene probes are in the columns and 'treatments' are in the rows, in a count matrix
WGCNA_matrix <- subCountData[,c(4,14:108)]

WGCNA_matrix <- WGCNA_matrix %>% 
  filter(., TargetName != "NegProbe-WTX")

rownames(WGCNA_matrix) <- WGCNA_matrix[, 1]

WGCNA_matrix <- WGCNA_matrix[, 2:96]

WGCNA_matrix <- as.matrix(WGCNA_matrix)

dim(WGCNA_matrix)

WGCNA_matrix_CA1_neun <- WGCNA_matrix[, grepl("neun", colnames(WGCNA_matrix)) & !grepl("rest", colnames(WGCNA_matrix)) & grepl("CA1", colnames(WGCNA_matrix)) | colnames(WGCNA_matrix) == "TargetName"]


WGCNA_matrix_CA1_rest <- WGCNA_matrix[, grepl("rest", colnames(WGCNA_matrix)) & !grepl("neun", colnames(WGCNA_matrix)) & grepl("CA1", colnames(WGCNA_matrix)) | colnames(WGCNA_matrix) == "TargetName"]


WGCNA_matrix_EC_neun <- WGCNA_matrix[, grepl("neun", colnames(WGCNA_matrix)) & !grepl("rest", colnames(WGCNA_matrix)) & grepl("EC", colnames(WGCNA_matrix)) | colnames(WGCNA_matrix) == "TargetName"]


WGCNA_matrix_EC_rest <- WGCNA_matrix[, grepl("rest", colnames(WGCNA_matrix)) & !grepl("neun", colnames(WGCNA_matrix)) & grepl("EC", colnames(WGCNA_matrix)) | colnames(WGCNA_matrix) == "TargetName"]


# Normalisation of matrix via cpm
#WGCNA_cpm <- edgeR::cpm(WGCNA_matrix)
#WGCNA_cpm_log <- edgeR::cpm(WGCNA_matrix, log =T)


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

#meta_df$age <- pathology$Month[pathology$Histology.no. %in% c("16/085", "16/095", "16/093", "16/083", "16/082", "16/092", "16/081", "16/099")]


rownames(meta_df) <- meta_df[,1]

# sorting for CA1 neun
meta_df <- meta_df %>%
  filter(., population == "neun" & grepl("CA1", rownames(meta_df)))
  

dds <- DESeqDataSetFromMatrix(round(WGCNA_matrix),
                              meta_df,
                              design = ~Group)
dds <- DESeq(dds)
# DEseq pipeline ends



#vsd_cpm <- varianceStabilizingTransformation(WGCNA_cpm)
# values are not integers

vsd <- varianceStabilizingTransformation(dds)
# what is this doing in the script?
# linear regression fit was substituted for a local regression fit by the function



wpn_vsd <- getVarianceStabilizedData(dds)

#wpn_vsd_cpm <- getVarianceStabilizedData(WGCNA_cpm)
#unable to find an inherited method for function ‘dispersionFunction’ for signature ‘object = "matrix"’

rv_wpn <- rowVars(wpn_vsd) 
#rv_cpm <- rowVars(WGCNA_cpm)
summary(rv_wpn)

q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
q5_wpn <- quantile( rowVars(wpn_vsd), .5) 
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= 95 quantile reduces dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
dim(expr_normalized)
# normalised expression obj showing 998 genes with 95q, and 4990 with 75q
# sometimes shows 1009??

#for cpm
#q75_cpm <- quantile( rowVars(WGCNA_cpm), .75)  # <= original
#q5_cpm <- quantile( rowVars(WGCNA_cpm), .5) 
#q95_cpm <- quantile( rowVars(WGCNA_cpm), .95)  # <= 95 quantile reduces dataset
#expr_normalized_cpm <- WGCNA_cpm[ rv_cpm > q75_cpm, ]
#dim(expr_normalized_cpm)
# same dim as dds above

#expr_normalized_cpm_log <- log(expr_normalized_cpm)

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
    title = "Normalized and 95 quantile Expression",
    x = "treatment",
    y = "normalized expression"
  )

# for cpm
#expr_normalized_cpm_log_df <- data.frame(expr_normalized_cpm_log) %>%
  mutate(
    Gene_id = row.names(expr_normalized_cpm_log)
  ) %>%
  pivot_longer(-Gene_id)

#expr_normalized_cpm_log_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "cpm Normalized",
    x = "treatment",
    y = "cpm normalized expression"
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
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.80, col = "red")

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

# trying power = 8




picked_power <- 8

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
#>            gene_id    colors
#> 1 AC149818.2_FG001      blue
#> 2 AC149829.2_FG003      blue
#> 3 AC182617.3_FG001      blue
#> 4 AC186512.3_FG001 turquoise
#> 5 AC186512.3_FG007 turquoise

write_delim(module_df,
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
MEs0 <- cbind(MEs0, meta_df)
#MEs0 <- left_join(MEs0, meta_df, by = "treatment")

# tidy & plot data
mME = MEs0 %>%
  select(- c("SegmentDisplayName", "Histology.no.", "Group", "grossRegion", "treatment")) %>%
  pivot_longer(cols = -population) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )





### Calculate correlation between traits and modules

MEs0 <- MEs0[,-which(colnames(MEs0)=="MEgrey")]
traits <- meta_df[,3:5]

#encoding traits
traits_numeric <- as.data.frame(lapply(traits, function(x) {
  if (is.factor(x) || is.character(x)) {
    as.numeric(as.factor(x))  # Converts factor levels to numbers
  } else {
    x  # Keep numeric columns as they are
  }
}))


moduleTraitCor = stats::cor(MEs0[,1:4],traits_numeric , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(MEs0))


png(paste0(path, title), width = 3000, height = 2000, res = 300)  # Increase resolution and dimensions
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

# Adjust margins and font sizes
par(mar = c(10, 10, 5, 5), cex.main = 1.5, cex.axis = 1, cex.lab = 0.7)  # Decrease margins and adjust font sizes

# Plot the heatmap

#pheatmap(moduleTraitCor, col = colours, main = "module trait correlation")


labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traits_numeric),
  yLabels = names(MEs0[1:4]),
  ySymbols = names(MEs0[1:4]),
  colorLabels = FALSE,
  colors = colorRampPalette(c("blue", "white", "red"))(50),  # Change color scheme if needed
  textMatrix = textMatrix,
  setStdMargins = TRUE,
  plotLegend = TRUE,
  cex.text = 0.7,  # Reduce label size
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)






MEs0 %>%
  select(- c("SegmentDisplayName", "Histology.no.", "Group", "grossRegion", "treatment")) %>%
  pivot_longer(cols = -population) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )




mME %>% ggplot(., aes(x= population, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")




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



# pick out a few modules of interest here
modules_of_interest = c("grey", "yellow")

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
  subset(colors %in% modules_of_interest)

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
            file = "edgelist.tsv",
            delim = "\t")

