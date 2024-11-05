
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
                               scale.fac=colData(spe_sub)$lib_size)
# changed scaling by nuclei count to lib size
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

