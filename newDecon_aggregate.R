BiocManager::install("SpatialDecon")
library(SpatialDecon)
library(tidyr)
library(dplyr)
library(reshape2)

# CPM normalise your spe object
spe_cpm <- geomxNorm(spe, method = "CPM")
spe_cpm <- addPerROIQC(spe_cpm, rm_genes = FALSE) # Add ROI QC

#Remove duplicates (this may differ per dataset)
duplicates <- c("D830030K20Rik", "Gm10406", "LOC118568634")
NormCountData <- assay(spe_cpm)
# NormCountData<- NormCountData[-which(NormCountData$TargetName %in% duplicates),]
# rownames(NormCountData) <- subCountData$ProbeDisplayName

norm <- NormCountData

negSet<- subCountData[grep("NegProbe",subCountData$TargetName),14:108]
negNames <- subCountData[grep("NegProbe",subCountData$TargetName),"ProbeDisplayName"]
rownames(negSet) <- negNames
norm <- rbind(norm,negSet)

# row.names(norm)[(nrow(norm) - 209):nrow(norm)] <- paste("NegProbe-WX",1:210)
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
colData(spe_cpm)[cellSubPropLong$variable,c("population","Group","grossRegion")]
cellSubPropLong <- cbind(cellSubPropLong,colData(spe_cpm)[cellSubPropLong$variable,c("population","Group","grossRegion")])
cellSubPropLong <- cellSubPropLong[order(cellSubPropLong$Group,cellSubPropLong$grossRegion),]
cellSubPropLong$variable <- as.character(cellSubPropLong$variable)
cellSubPropLong$variable <- factor(cellSubPropLong$variable, levels = unique(cellSubPropLong$variable))
cellSubPropLong$Group <- factor(cellSubPropLong$Group,levels = c("WW","CC"))


pdf("NewAllenDecon.pdf",width = 6, height = 8)
ggplot(cellSubPropLong,aes(y = `sum(value)`, x = as.factor(BroadCell), color = as.factor(Group)))+
  geom_boxplot()+ 
  geom_point(position=position_jitterdodge(),alpha = 0.7,size = 0.6)+
  facet_grid(population~grossRegion)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0.01))+
  theme_bw()+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(color = "Transgene Group",x = "Broad Cell Type (Allen Brain Atlas)", y = "Proportion")+
  scale_color_brewer(palette = "Set2",labels = c("Wildtype","Transgenic"))
dev.off() 
