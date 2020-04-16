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

WT <- readMM('TREM2/Data/matrix_WT.mtx')
rownames(WT) <- readLines('TREM2/Data/geneNames_WT.txt')
colnames(WT) <- readLines('TREM2/Data/barcodes_WT.txt')
WT <- scQC(WT)
WT <- WT[rowMeans(WT != 0) > 0.05,]
WT <- WT[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(WT), ignore.case = TRUE),]
WT <- scTenifoldKnk(WT, gKO = 'Trem2')
save(WT, file = 'TREM2_GSE130626.RData')

