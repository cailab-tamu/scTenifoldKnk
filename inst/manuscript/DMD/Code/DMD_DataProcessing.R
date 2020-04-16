setwd('/data/dcosorioh/manuscript/')

library(Matrix)
scTenifoldKnk <- function(countMatrix, gKO = NULL){
  set.seed(1)
  WT <- scTenifoldNet::makeNetworks(countMatrix, q = 0.9)
  set.seed(1)
  WT <- scTenifoldNet::tensorDecomposition(WT)
  WT <- as.matrix(WT$X)
  #KO <- rCUR::CUR(WT, sv = RSpectra::svds(WT, 5))
  #C <- KO@C
  #C[gKO,] <- 0
  #KO <- C %*% KO@U %*% KO@R
  KO <- WT
  KO[gKO,] <- 0
  set.seed(1)
  MA <- scTenifoldNet::manifoldAlignment(WT, KO)
  set.seed(1)
  DR <- scTenifoldNet::dRegulation(MA)
  outputList <- list()
  outputList$WT <- WT
  outputList$KO <- KO
  outputList$diffRegulation <- DR
  return(outputList)
}
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')

GSM4116571 <- readMM('DMD/Data/GSM4116571_Qc_matrix.mtx.gz')
rownames(GSM4116571) <- read.csv('DMD/Data/GSM4116571_Qc_features.tsv.gz', sep = '\t', header = FALSE, stringsAsFactors = FALSE)[,2]
colnames(GSM4116571) <- readLines('DMD/Data/GSM4116571_Qc_barcodes.tsv.gz')
GSM4116571 <- scQC(GSM4116571)
GSM4116571 <- GSM4116571[rowMeans(GSM4116571 != 0) > 0.05,]
GSM4116571 <- GSM4116571[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(GSM4116571), ignore.case = TRUE),]
GSM4116571 <- scTenifoldKnk(GSM4116571, gKO = 'Dmd')
save(GSM4116571, file = 'DMD_GSM4116571.RData')

