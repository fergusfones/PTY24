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

metaData <- read.csv("preProcessedMetaData_filtered.csv",header = T)
countDataa <- read.csv("preProcessedCountData_filtered.csv",header = T,)
pathology <- read.csv("PathologyPG5.csv", header = T,)
rownames(pathology) <- pathology$Histology.no.
rownames(metaData) <- metaData$SegmentDisplayName
colnames(countDataa)[14:108] <- gsub("X","",colnames(countDataa)[14:108])

duplicates <- c("D830030K20Rik", "Gm10406", "LOC118568634")
subCountData <- countDataa[-which(countDataa$TargetName %in% duplicates),]  
rownames(subCountData) <- subCountData$ProbeDisplayName


# need to add covaraties before running WGCNA
# do i have to control for cell types and correct the data like we did in the GeoMX workflow?
# just use neun samples???


# need to transpose data, so gene probes are in the coloumns and 'treatments' are in the rows, in a count matrix
WGCNA_matrix <- subCountData[,c(4,14:108)]

WGCNA_matrix <- WGCNA_matrix %>% 
  filter(., TargetName != "NegProbe-WTX")

rownames(WGCNA_matrix) <- WGCNA_matrix[, 1]
WGCNA_matrix <- WGCNA_matrix[, 2:96]

WGCNA_matrix <- as.matrix(WGCNA_matrix)

# Normalising data with DESeq
meta_df <- read.csv("preProcessedMetaData_filtered.csv",header = T)
meta_df <- meta_df[, c(6, 36, 38 )]
rownames(meta_df) <- meta_df[,1]

dds <- DESeqDataSetFromMatrix(round(WGCNA_matrix),
                              meta_df,
                              design = ~Group + population)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds)
# linear regression fit was substituted for a local regression fit by the function

wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)

q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
dim(expr_normalized)
# normalised expression obj showing 998 genes now...

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

# starting WGCNA

input_mat <- t(expr_normalized)

# Allowing multi-threading
allowWGCNAThreads() 

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  expr_normalized,             # <= Input data
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
abline(h = 0.90, col = "red")
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
