# The code performs a WGCNA analysis
# 
# Made using this tutorial : https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0
# And with help from Giulia Pegoraro, University of Exeter Medical School

###########################################
setwd("/Users/fergusfones/Desktop/Nanostring/")
getwd()


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
#install.packages("RColorBrewer")
library(RColorBrewer)
library(magrittr)
library(clusterProfiler)
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(purrr)
BiocManager::install("VennDiagram")
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



# do i have to control for covariates like cell types and correct the data like we did in the GeoMX workflow?


# need to transpose data, so gene probes are in the columns and 'treatments' are in the rows, in a count matrix
WGCNA_matrix <- subCountData[,c(4,14:108)]

WGCNA_matrix <- WGCNA_matrix %>% 
  filter(., TargetName != "NegProbe-WTX")

rownames(WGCNA_matrix) <- WGCNA_matrix[, 1]

WGCNA_matrix <- WGCNA_matrix[, 2:96]

WGCNA_matrix <- as.matrix(WGCNA_matrix)

dim(WGCNA_matrix)
# 19959, 95

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


dds <- DESeqDataSetFromMatrix(round(WGCNA_matrix),
                              meta_df,
                              design = ~Group)
dds <- DESeq(dds)

# DEseq pipeline ends
save(dds, file = "/Users/fergusfones/Desktop/Nanostring/dds.Rdata")
load(file = "/Users/fergusfones/Desktop/Nanostring/dds.Rdata")

vsd <- varianceStabilizingTransformation(dds)
# retains metadata, can be used for PCA, visualisation, further processing
# linear regression fit was substituted for a local regression fit by the function



wpn_vsd <- getVarianceStabilizedData(dds)


rv_wpn <- rowVars(wpn_vsd) 
summary(rv_wpn)

q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
q5_wpn <- quantile( rowVars(wpn_vsd), .5) 
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= 95 quantile reduces dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
dim(expr_normalized)
# whole normalised expression obj showing 998 genes with 95q, and 4990 with 75q



#
#
#
#
#
#
#
#
# Subset after normalisation
#
#
#
#
#
#
#

expr_normalized_CA1_neun <- expr_normalized[, grepl("neun", colnames(expr_normalized)) & !grepl("rest", colnames(expr_normalized)) & grepl("CA1", colnames(expr_normalized)) | colnames(expr_normalized) == "TargetName"]


expr_normalized_CA1_rest <- expr_normalized[, grepl("rest", colnames(expr_normalized)) & !grepl("neun", colnames(expr_normalized)) & grepl("CA1", colnames(expr_normalized)) | colnames(expr_normalized) == "TargetName"]


expr_normalized_EC_neun <- expr_normalized[, grepl("neun", colnames(expr_normalized)) & !grepl("rest", colnames(expr_normalized)) & grepl("EC", colnames(expr_normalized)) | colnames(expr_normalized) == "TargetName"]


expr_normalized_EC_rest <- expr_normalized[, grepl("rest", colnames(expr_normalized)) & !grepl("neun", colnames(expr_normalized)) & grepl("EC", colnames(expr_normalized)) | colnames(expr_normalized) == "TargetName"]



save(expr_normalized_CA1_neun, file = "/Users/fergusfones/Desktop/Nanostring/expr_normalized_CA1_neun.Rdata")
load(file = "/Users/fergusfones/Desktop/Nanostring/expr_normalized_CA1_neun.Rdata")

save(expr_normalized_CA1_rest, file = "/Users/fergusfones/Desktop/Nanostring/expr_normalized_CA1_rest.Rdata")
load(file = "/Users/fergusfones/Desktop/Nanostring/expr_normalized_CA1_rest.Rdata")

save(expr_normalized_EC_neun, file = "/Users/fergusfones/Desktop/Nanostring/expr_normalized_EC_neun.Rdata")
load(file = "/Users/fergusfones/Desktop/Nanostring/expr_normalized_EC_neun.Rdata")

save(expr_normalized_EC_rest, file = "/Users/fergusfones/Desktop/Nanostring/expr_normalized_EC_rest.Rdata")
load(file = "/Users/fergusfones/Desktop/Nanostring/expr_normalized_EC_rest.Rdata")


expr_normalized_df <- data.frame(expr_normalized_CA1_neun) %>%
  mutate(
    Gene_id = row.names(expr_normalized_CA1_neun)
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
  
  
  cor_matrix <-  t(expr_normalized_EC_rest)
  cor_matrix <- cor(cor_matrix, method = "spearman")
  
  colours <- colorRampPalette(brewer.pal(9, "Blues"))(225)
  
  pheatmap(cor_matrix, col = colours, main = "gene expression correlation", show_colnames = F, show_rownames = F)
  
  
  
  # hierarchical clustering 
  
htree <- hclust(dist(t(expr_normalized_EC_rest)))
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

input_mat <- t(expr_normalized_CA1_neun)

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

# try using the estimate power function

picked_power <- 3

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
MEs0 <- cbind(MEs0, meta_df_CA1_neun)


### Calculate correlation between traits and modules

MEs0 <- MEs0[,-which(colnames(MEs0)=="MEgrey")]

traits <- meta_df_CA1_neun %>%
  select(- c("SegmentDisplayName", "Histology.no.", "population", "grossRegion"))

#encoding traits
traits_numeric <- ifelse(traits == "WW", 0 ,1) 

#traits_numeric <- traits %>%
#  mutate(
 # CC = ifelse(traits$Group == "CC", 1, 0),
#  WW = ifelse(traits$Group == "WW", 1, 0)
#)

#traits_numeric <- traits_numeric[, 2:3]

moduleTraitCor = stats::cor(MEs0[,1:4],traits_numeric , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(MEs0))


png(paste0(path, title), width = 3000, height = 2000, res = 300)  # Increase resolution and dimensions
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

# tidy & plot data
#mME = MEs0 %>%
  select(- c("SegmentDisplayName", "Histology.no.", "population", "grossRegion", "treatment")) %>%
  pivot_longer(cols = -Group) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )


#mME %>% ggplot(., aes(x= Group, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships_CA1_rest", y = "Modules", fill="corr")










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
modules_of_interest = c("green", "yellow", "black", "turquoise")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
expr_normalized_CA1_neun[1:5,1:10]

subexpr = expr_normalized_CA1_neun[submod$gene_id,]

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



######
#
# The function performs KEGG enrichment on significant results from two-way ANOVA
# module_df: dataframe of modules resulting from WGCNA
# path: character containing path in which to save the images
#
path <- ("/Users/fergusfones/Desktop/Nanostring/KEGG/")
######

KEGG_module <- function(module_df, path){
  col_enrich_kegg <- lapply(unique(module_df$colors), function(color){
    #gene_list <- strsplit(module_df[which(module_df$colors == color),"gene_id"],"\\.")
    #gene_list <- lapply(1:length(gene_list), function(x){ gene_list[[x]][[1]]})
    gene_list <-  module_df[which(module_df$colors == color),"gene_id"]
    gene.df <- bitr(gene_list, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Mm.eg.db)
    kegen <- enrichKEGG(gene     = gene.df$ENTREZID,
                        organism     = "mmu",
                        pvalueCutoff = 0.05)
  })
  names(col_enrich_kegg) <- unique(module_df$colors)
  col_enrich_kegg <- col_enrich_kegg %>% keep( ~ nrow(.) !=0 )
  lapply(names(col_enrich_kegg), function(i){
    file_path <- paste0(path,i,"_enrich.png")
    barplot((col_enrich_kegg[[i]]), showCategory = 10)
    ggsave(file_path)
  })
}

results <- KEGG_module(module_df, path)

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
# dont think this is working properly

expr_of_interest = expr_normalized_CA1_neun[genes_of_interest$gene_id,]
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
            file = "/Users/fergusfones/Desktop/Nanostring/edge_list.tsv",
            delim = "\t")

#
#
#
#
#
venn.diagram(x = list(module_df$gene_id[module_df$colors == "turquoise"], rownames(de_genes_toptable_CA1)),
             category.names = c("WGCNA turq module", "DGE genes"),
             filename = '#14_venn_diagramm.png',
             output=TRUE
)

WGCNA_turq <-module_df$gene_id[module_df$colors == "turquoise"]
toptable_CA1<- rownames(de_genes_toptable_CA1)

write_lines(WGCNA_turq, file = "/Users/fergusfones/Desktop/Nanostring/WGCNA_turq.txt")
write_lines(toptable_CA1, file = "/Users/fergusfones/Desktop/Nanostring/toptable_CA1.txt")

turq_wgcna<- intersect(WGCNA_turq,toptable_CA1)





#enrichment

# turqouise module of CA1_nuen contains Gfap


library(enrichR)

dbs <- listEnrichrDbs()

head(dbs)


dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023")

turq_enriched <- enrichr(WGCNA_turq, dbs)

plotEnrich(turq_enriched[[1]], showTerms = 20, numChar = 120, y = "Count", orderBy = "P.value")
