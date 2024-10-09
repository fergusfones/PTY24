# This function takes a spe object and a pre-calculated set of principle components and tests/plots
# the results of correlation testing for each individual PC. It requires the user to have run the
# scater::runPCA() and the reducedDim() functions already

# Package requirements

corPCtest <- function(speOb, pcaOb, nPC, colLabs){

#speOb = spe object (with pcs calculated)
#pcaOb = principle component matrix from the reducedDim function
#nPC = number of principle component
#colLabs = vector of column names in the colData from the spe data to test  
  
    
require(ggplot2)
require(tidyverse)
require(dplyr)
require(fastDummies)  


# Check data objects
if(unique(rownames(colData(speOb)) == rownames(pcaOb)) == TRUE){
  print("Data matches between spe colData and PCA data")
}else{
  print("ERROR: Mismatched spe colData and PCA data, check data is consistent")
}


tempPheno <- colData(speOb)[,colLabs]
rmNames <- c()
for(dat in colnames(tempPheno)){
  if(class(tempPheno[,dat]) == "character"){
    tempPheno[,dat] <- as.factor(tempPheno[,dat])
    if(length(levels(tempPheno[,dat])) > 2){
      print(paste(dat,"requires dummy variable transformation"))
      tempPheno <- cbind(tempPheno,fastDummies::dummy_cols(tempPheno[,dat])[,-1])
      tempPheno <- tempPheno[,-which(colnames(tempPheno) == dat)]
    }else{
      tempPheno[,dat] <- as.numeric(tempPheno[,dat])
    }
  }
}

tempPC <- pcaOb[,1:nPC]

# Generate correlation matrices and test significance
corTes <- cor(as.matrix(tempPheno),as.matrix(tempPC))

corTes <- melt(corTes)
corTes$P <- NA
for(tes in 1:nrow(corTes)){
  phenoVar <- corTes[tes,"Var1"]
  PCvar <- corTes[tes,"Var2"]
  corTes$P[tes] <- cor.test(tempPheno[,phenoVar],tempPC[,PCvar])$p.value
}

corTes$Plab <- NA
corTes[which(corTes$P < 0.05),"Plab"] <- "*"
corTes[which(corTes$P < 0.01),"Plab"] <- "**"
corTes[which(p.adjust(corTes$P) < 0.05),"Plab"] <- "***"

corTes$lab <- paste(signif(corTes$value,2), "\n",corTes$Plab)
corTes[which(is.na(corTes$Plab)),"lab"] <- NA

corHeat <- ggplot(corTes,aes(x = Var1,y = Var2, fill = abs(value)))+
  geom_tile()+
  scale_fill_viridis_c()+
  geom_text(aes(label = lab))+
  xlab("Variable")+
  ylab("PC")+
  labs(fill = "Absolute Correlation")+
  scale_y_discrete(limits = rev(levels(corTes$Var2)))
  
print(corHeat)
corTes <- corTes[,c(1:4)]
colnames(corTes) <- c("Variable","PC","Correlation","P-value")
return(corTes)
}


