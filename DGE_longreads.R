long_reads<- read.delim("/Users/fergusfones/Desktop/Nanostring/rTg4510SCN_collapsed_RulesFilter_result_classification_counts.txt")
phenotypes <- read.csv("/Users/fergusfones/Desktop/Nanostring/SCNPhenotype.csv")

save(long_reads, file = "/Users/fergusfones/Desktop/Nanostring/rTg4510SCN_collapsed_RulesFilter_result_classification_counts.txt")
load("/Users/fergusfones/Desktop/Nanostring/rTg4510SCN_collapsed_RulesFilter_result_classification_counts.txt")

save(phenotypes, file = "/Users/fergusfones/Desktop/Nanostring/SCNPhenotype.csv")
load("/Users/fergusfones/Desktop/Nanostring/SCNPhenotype.csv")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

BiocManager::install("Gmisc")

library(DESeq2)
library(tibble)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(standR)
library(Gmisc, quietly = TRUE)
library(grid)
library(glue)
library(dplyr)
library(edgeR)

identical(long_reads$NeuN65_mapped, long_reads$NeuN65._mapped)

long_reads <- long_reads[ -69]

# sum the counts of transcript by genes
cts <- long_reads %>%
  group_by(associated_gene) %>%
  summarise((across(58:72, sum, na.rm = TRUE)))

(unique(long_reads$associated_gene))
# looks like the genes have been summed right

#put the genes into the rows

cts <- tibble::column_to_rownames(cts, var = "associated_gene")

# subset Neun
cts <- cts[, c( "NeuN15_mapped", "NeuN19_mapped" ,"NeuN3_mapped" ,"NeuN65_mapped" ,"NeuN69_mapped" ,"NeuN72_mapped" ,"NeuN77_mapped", "NeuN7_mapped")]


#subset Neun in colData
phenotypes <- phenotypes[9:16,]

# making the DEseqDataset object, i think its like a summarised experiment object...
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = phenotypes,
                              design= ~group)

keep <- rowSums(counts(dds)) >= 20
# upped the threshold to 20

dds <- dds[keep,]


DE_dds <- DESeq(dds)
#Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : every gene contains at least one zero, cannot compute log geometric means
# could add 1 to all counts, however this will reduce effect size of counts that are small...

cts_normalised <- edgeR::cpm(cts)
cts_normalised_log <- edgeR::cpm(cts, log = TRUE)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = phenotypes,
                              design= ~group)





# have no idea what i was trying here
##p <- cbind(reshape2::melt(sizeFactors(dds)), reshape2::melt(colSums(counts(dds)))) %>% 
  `colnames<-`(c("sizefactors", "nreads")) %>% 
  tibble::rownames_to_column("sample") %>% 
  mutate(sample = str_remove(sample,"ont_")) %>%
  ggplot(., aes(x = sizefactors, y = nreads)) + geom_point() +
  geom_label_repel(aes(label = sample), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  theme_bw() + labs(y = "Number of reads", x = "Size Factors")
p  


# rlog is the recommended method for experiments with n<20
# normalization to stabilize variance (regularized logarithm)
rld <- rlog(dds, blind = FALSE)





ggplot(dds,
       aes(x = group, y =counts(dds) , fill = group) +
         geom_bar(stat = "identity", position = "dodge") +
         labs(
         x = "Age and Group",
         y = "Sample Number (n)",
         title = "Overview of the Dataframe") +
         theme_minimal() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
         scale_fill_manual(values = c("Tg" = "skyblue", "WT" = "orange")))






# trying to make a boxplot

org_cohort <- boxGrob(glue("mice",
                           "n = 8",
                           .sep = "\n"))

included <- boxGrob(glue("NeuN+ cells",
                         "n = 8",
                         .sep = "\n"))

excluded <- boxGrob(glue("Excluded DN cells",
                         "n = 8",
                         .sep = "\n"),
                    just = "left")

two_months <- boxGrob(glue("2 months",
                      "n = 4",
                      .sep = "\n"))

eight_months <- boxGrob(glue("8 months",
                      "n = 4",
                      .sep = "\n"))

WT_group <- boxGrob(glue("WT",
                         "n = 2",
                         .sep = "\n"))

Tg_group <- boxGrob(glue("Tg",
                         "n = 2",
                         .sep = "\n"))

grid.newpage()
vert <- spreadVertical(org_cohort,
                       included = included,
                       grps = two_months)
grps <- alignVertical(reference = vert$grps,
                      two_months, eight_months) %>%
  spreadHorizontal()
vert$grps <- NULL

excluded <- moveBox(excluded,
                    x = .1,
                    y = coords(vert$included)$top + distance(vert$org_cohort, vert$gprs, half = TRUE, center = FALSE))

for (i in 1:(length(vert) - 1)) {
  connectGrob(vert[[i]], vert[[i + 1]], type = "vert") %>%
    print
}
connectGrob(vert$included, grps[[1]], type = "N")
connectGrob(vert$included, grps[[2]], type = "N")

connectGrob(vert$eligible, excluded, type = "L")

# Print boxes
vert
grps
excluded
       
##################

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




Pfkm <- long_reads[long_reads$associated_gene == "Pfkm",] %>%
  select(contains("mapped")) %>%
  summarise(across(1:16, ~ sum(.x, na.rm = TRUE)))

Pfkm %>% reshape2::melt(variable.name = "sample",value.name = "reads") %>% 
  mutate(sample_id = stringr::word(sample,c(1), sep = fixed("_"))) %>% 
  merge(., phenotypes, by = "sample") %>% 
  mutate(group = factor(group, levels = c("WT","TG"))) %>%
  ggplot(., aes(x = group, y = reads, colour = cell)) + geom_point(size = 3) +
  facet_grid(~age)




Cdkn2aip <- long_reads[long_reads$associated_gene == "Cdkn2aip",] %>%
  select(contains("mapped")) %>%
  summarise(across(1:16, ~ sum(.x, na.rm = TRUE)))

Cdkn2aip %>% reshape2::melt(variable.name = "sample",value.name = "reads") %>% 
  mutate(sample_id = stringr::word(sample,c(1), sep = fixed("_"))) %>% 
  merge(., phenotypes, by = "sample") %>% 
  mutate(group = factor(group, levels = c("WT","TG"))) %>%
  ggplot(., aes(x = group, y = reads, colour = cell)) + geom_point(size = 3) +
  facet_grid(~age)



Clec2j <- long_reads[long_reads$associated_gene == "Cdkn2aip",] %>%
  select(contains("mapped")) %>%
  summarise(across(1:16, ~ sum(.x, na.rm = TRUE)))

Clec2j %>% reshape2::melt(variable.name = "sample",value.name = "reads") %>% 
  mutate(sample_id = stringr::word(sample,c(1), sep = fixed("_"))) %>% 
  merge(., phenotypes, by = "sample") %>% 
  mutate(group = factor(group, levels = c("WT","TG"))) %>%
  ggplot(., aes(x = group, y = reads, colour = cell)) + geom_point(size = 3) +
  facet_grid(~age)

