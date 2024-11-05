# substetting before RUV corrections

spe_CA1_neun <- spe[,colData(spe)$grossRegion == "CA1" &
                              colData(spe)$population == "neun"]

spe_CA1_rest <- spe[,colData(spe)$grossRegion == "CA1" &
                              colData(spe)$population == "rest"]

spe_EC_neun <- spe[,colData(spe)$grossRegion == "EC" &
                             colData(spe)$population == "neun"]

spe_EC_rest <- spe[,colData(spe)$grossRegion == "EC" &
                             colData(spe)$population == "rest"]

spe_CA1_neun <- findNCGs(spe_CA1_neun, batch_name = "SlideName", top_n = 300)
spe_CA1_rest <- findNCGs(spe_CA1_rest, batch_name = "SlideName", top_n = 300)
spe_EC_neun <- findNCGs(spe_EC_neun, batch_name = "SlideName", top_n = 300)
spe_EC_rest <- findNCGs(spe_EC_rest, batch_name = "SlideName", top_n = 300)

for(i in seq(5)){
  spe_CA1_neun_ruv_post <- geomxBatchCorrection(spe_CA1_neun, factors = c("Group","grossRegion"), 
                                  NCGs = metadata(spe_CA1_neun)$NCGs, k = i)
  
  print(plotPairPCA(spe_CA1_neun_ruv_post, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
}

for(i in seq(5)){
  spe_CA1_rest_ruv_post <- geomxBatchCorrection(spe_CA1_rest, factors = c("Group","grossRegion"), 
                                                NCGs = metadata(spe_CA1_rest)$NCGs, k = i)
  
  print(plotPairPCA(spe_CA1_rest_ruv_post, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
}

for(i in seq(5)){
  spe_EC_neun_ruv_post <- geomxBatchCorrection(spe_EC_neun, factors = c("Group","grossRegion"), 
                                                NCGs = metadata(spe_EC_neun)$NCGs, k = i)
  
  print(plotPairPCA(spe_EC_neun_ruv_post, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
}

for(i in seq(5)){
  spe_EC_rest_ruv_post <- geomxBatchCorrection(spe_EC_rest, factors = c("Group","grossRegion"), 
                                                NCGs = metadata(spe_EC_rest)$NCGs, k = i)
  
  print(plotPairPCA(spe_EC_rest_ruv_post, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
}


plotPairPCA(spe_CA1_neun_ruv_post, assay = 2, color = Group, title = "spe_CA1_neun_ruv_post", n_dimension = 4)
plotPairPCA(spe_CA1_rest_ruv_post, assay = 2, color = Group, title = "spe_CA1_rest_ruv_post", n_dimension = 4)
plotPairPCA(spe_EC_neun_ruv_post, assay = 2, color = Group, title = "spe_EC_neun_ruv_post", n_dimension = 4)
plotPairPCA(spe_EC_rest_ruv_post, assay = 2, color = Group, title = "spe_EC_rest_ruv_post", n_dimension = 4)




