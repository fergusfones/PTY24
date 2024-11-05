parameters <- c("grossRegion", "SlideName", "Group", "population")
list = list(spe$Group, spe$SlideName, spe$population, spe$grossRegion)
par(mfrow = c(2, 2))

loop.vector <- 1:4

for (i in seq(1:4)){
  print(i)
  png(paste0(i,".png"))
  plotPairPCA(spe[,i], precomputed = pca_results, n_dimension = 4)
  dev.off()
  }

lapply(parameters, function(x){
  drawPCA(spe, precomputed = pca_results, color = x)
})


# try and alter this
for(i in seq(5)){
  spe_ruv <- geomxBatchCorrection(spe, factors = "Type", 
                                  NCGs = metadata(spe)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Type, title = paste0("k = ", i)))
  
}
