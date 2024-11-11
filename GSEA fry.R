## GSEA

library(msigdb)
library(GSEABase)

msigdb_mm <- getMsigdb(version = '7.2', org ='mm', id='SYM')
msigdb_mm <- appendKEGG(msigdb_mm)

sc <- listSubCollections(msigdb_mm)

gsc <- c(subsetCollection(msigdb_mm, 'c2', sc[grepl("^CP:",sc)]),
         subsetCollection(msigdb_mm, 'c5', sc[grepl("^GO:",sc)])) %>%
  GeneSetCollection()


fry_indices <- ids2indices(lapply(gsc, geneIds), rownames(v), remove.empty = FALSE)
names(fry_indices) <- sapply(gsc, setName)

gsc_category <- sapply(gsc, function(x) bcCategory(collectionType(x)))
gsc_category <- gsc_category[sapply(fry_indices, length) > 5]

gsc_subcategory <- sapply(gsc, function(x) bcSubCategory(collectionType(x)))
gsc_subcategory <- gsc_subcategory[sapply(fry_indices, length) > 5]

fry_indices <- fry_indices[sapply(fry_indices, length) > 5]

names(gsc_category) = names(gsc_subcategory) = names(fry_indices)

fry_indices_cat <- split(fry_indices, gsc_category[names(fry_indices)])
fry_res_out <- lapply(fry_indices_cat, function (x) {
  limma::fry(v, index = x, design = design, contrast = contr.matrix[,1], robust = TRUE)
})

post_fry_format <- function(fry_output, gsc_category, gsc_subcategory){
  names(fry_output) <- NULL
  fry_output <- do.call(rbind, fry_output)
  fry_output$GenesetName <- rownames(fry_output)
  fry_output$GenesetCat <- gsc_category[rownames(fry_output)]
  fry_output$GenesetSubCat <- gsc_subcategory[rownames(fry_output)]
  return(fry_output)
}

fry_res_sig <- post_fry_format(fry_res_out, gsc_category, gsc_subcategory) %>%
  as.data.frame()
#  filter(FDR < 0.05) 

library(ggplot2)
library(dplyr)

fry_res_sig %>%
  arrange(FDR) %>%
  filter(Direction == "Up") %>%
  .[seq(20),] %>%
  mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
  ggplot(aes(GenesetName, -log(FDR))) +
  geom_bar(stat = "identity", fill = "red") +
  theme_bw() +
  coord_flip() +
  ggtitle("Up-regulated")

fry_res_sig %>%
  arrange(FDR) %>%
  filter(Direction == "Down") %>%
  .[seq(20),] %>%
  mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
  ggplot(aes(GenesetName, -log(FDR))) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_bw() +
  coord_flip() +
  ggtitle("Down-regulated")
