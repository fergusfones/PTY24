setwd("/Users/fergusfones/Desktop/Nanostring/")
getwd()

library(devtools)
install_github("wjawaid/enrichR")
library(enrichR)

dbs <- listEnrichrDbs()

head(dbs)


dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023")

### enriched subsets

rownames(de_results)
x <- rownames(top_n(de_results, 25))

spe_CA1_neun_ruv_post_enrich <- enrichr(c("Mapt",     "Car4"   ,       "Slc35f1"   ,    "Tbr1"   ,       "Hrk"         ,  "Lmbrd2"     ,   "Gfap"    ,      "Crym"     ,     "Cap1"  ,        "Celf2"     ,    "Ptprn2"  ,      "Tspan5"     ,   "App"   ,        "Cyfip2"   ,     "Acap2" ,"Grina"     ,    "Tuba1a"    ,    "Txndc15" ,      "Zbtb18"      ,  "Atp2b2"   ,     "Sncb"     ,     "Epha5"      ,   "Pcsk2"     ,    "Marcksl1"   ,   "Stim2"   ,      "Igfbp4"   ,     "Tmsb4x"    ,    "Rps24"       ,  "Slc30a3" ,      "Wasf1" ), dbs)
spe_CA1_rest_ruv_post_enrich <- enrichr(c("Mapt"  ,    "Gfap"   ,       "Ctsd"      ,    "C4b"       ,    "Ubb"    ,       "Cyfip2"    ,    "Nrgn"    ,      "Ctss"     ,     "Itm2b"     ,    "Serpina3n"  ,   "Hrk"     ,      "Ncdn"     ,     "Calm2"   , "Zbtb18"  ,      "C1qa"      ,    "Ppp3ca"    ,    "Mt1"     ,      "Fkbp1a"   ,     "Scn1b"    ,     "Cfl1"       ,   "App"    ,       "Hexb"       ,   "Ppp3r1"   ,     "B2m"      ,     "Clu"     ,      "Capza2"   ,"Actb"    ,      "Thy1"      ,    "Actr2"   ,      "Mt2"), dbs)
spe_EC_neun_ruv_post_enrich <- enrichr(c("Mapt"     ,     "Gfap"      ,    "Ptprn2"    ,    "Cst3"      ,    "Apoe"      ,    "Aldoa"     ,    "Smc3"      ,    "Ptgds"    ,     "Fer1l6"      ,  "Fgf14"  ,       "Esyt2"  ,      "Glul"  ,        "Slc17a7"  , "Ndrg2"    ,     "C1qa"  ,        "Aldoc"     ,    "Sncb"    ,      "Cbx8"      ,    "Mdh1"      ,    "Zwint"     ,    "Syn2"    ,      "Nap1l5"     ,   "Atp1a2"    ,    "Mt1"    ,       "Gapdh"     ,    "Nsf" , "Tmc5"    ,      "Stmn1"    ,     "Rab3a"   ,      "Atp6v0c"), dbs)
spe_EC_rest_ruv_post_enrich <- enrichr(c("Mapt"   ,     "Gfap"    ,      "Ctsd"     ,     "Ptprn2"      ,  "Camk2n1"    ,   "C1qa"    ,      "C4b"    ,       "Ctss"      ,    "Id3"   ,        "Fer1l6"    ,    "Nrsn1"     ,    "Cst3"   ,       "Id2"    ,"Ptgds"   ,      "Ighm"    ,      "Smc3"     ,     "Cplx1"      ,   "Ankrd12"   ,    "Cbx8"    ,      "Tmsb4x"  ,      "Plp1"     ,     "C1qc"    ,      "Clu"      ,     "Cox8a"    ,     "Gm52800"     ,  "Nfasc"   , "Cd81"   ,       "Cd9"      ,     "Tyrobp"    ,    "C1qb"), dbs)

plotEnrich(spe_CA1_neun_ruv_post_enrich[[1]], showTerms = 20, numChar = 120, y = "Count", orderBy = "P.value")
plotEnrich(spe_CA1_rest_ruv_post_enrich[[3]], showTerms = 20, numChar = 120, y = "Count", orderBy = "P.value")
plotEnrich(spe_EC_neun_ruv_post_enrich[[3]], showTerms = 20, numChar = 120, y = "Count", orderBy = "P.value")
plotEnrich(spe_EC_rest_ruv_post_enrich[[3]], showTerms = 20, numChar = 120, y = "Count", orderBy = "P.value")



#whole gene sets

install.packages("tibble")
library(tibble)
library(tidyverse)

de_genes_toptable <- topTable(fit2, coef = 1, sort.by = "p", n = Inf,p.value = 0.05,adjust.method = "fdr",lfc = 0.5)
offTargets <- c("Wdr60", "Esyt2", "Ncapg2","Ptprn2","Fgf14")

de_genes_toptable <- de_genes_toptable[-which(rownames(de_genes_toptable) %in% offTargets),]

whole_spe_CA1_neun <- enrichr(rownames(de_genes_toptable), dbs)

plotEnrich(whole_spe_CA1_neun[[3]], showTerms = 20, numChar = 120, y = "Count", orderBy = "P.value")

