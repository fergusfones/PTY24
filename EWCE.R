####### Script to test the enrichment of DEG signature in a mouse reference single cell
####### dataset using EWCE ######
getwd()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EWCE")

library(EWCE)
library(ewceData)

mouseCTD <- ewceData::ctd()

#Let's use an example topTable for one brain region
#subset significant DEG's from your toptables 
de_genes_toptable_EC_rest <- topTable(fit2_EC_rest, coef = 1, sort.by = "p", n = Inf,p.value = 0.05,adjust.method = "fdr",lfc = 0.5)
offTargets <- c("Wdr60", "Esyt2", "Ncapg2","Ptprn2","Fgf14")
#Manually remove any known off-targets
de_genes_toptable_EC_rest <- de_genes_toptable_EC_rest[-which(rownames(de_genes_toptable_EC_rest) %in% offTargets),]

#finally subset the list of names for your DEGs
topGenes <- rownames(de_genes_toptable_EC_rest)


# This tests for enrichment in the mouse single cell dataset
# Start out with 100 for quick running but when we run the analysis "finally" make this 10,000 reps
results <- EWCE::bootstrap_enrichment_test(sct_data = mouseCTD,
                                           sctSpecies = "mouse",
                                           genelistSpecies = "mouse",
                                           hits = topGenes, 
                                           reps = 100,
                                           annotLevel = 1)

saveRDS(results, file = "/Users/fergusfones/Desktop/Nanostring/results_CA1_rest.Rdata")
results_CA1_rest <- readRDS(file = "/Users/fergusfones/Desktop/Nanostring/results_CA1_rest.Rdata")

saveRDS(results, file = "/Users/fergusfones/Desktop/Nanostring/results_EC_rest.Rdata")
results_EC_rest <- readRDS(file = "/Users/fergusfones/Desktop/Nanostring/results_EC_rest.Rdata")

#Take out the enrichment results
ewceRes <- results_EC_rest$results


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
