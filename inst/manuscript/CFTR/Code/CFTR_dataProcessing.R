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

SRS3161261 <- readMM('CFTR/Data/matrix_SRS3161261.mtx')
rownames(SRS3161261) <- readLines('CFTR/Data/geneNames_SRS3161261.txt')
colnames(SRS3161261) <- readLines('CFTR/Data/barcodes_SRS3161261.txt')
SRS3161261 <- scQC(SRS3161261)
SRS3161261 <- SRS3161261[rowMeans(SRS3161261 != 0) > 0.05,]
SRS3161261 <- SRS3161261[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(SRS3161261), ignore.case = TRUE),]
SRS3161261 <- scTenifoldKnk(SRS3161261)
save(SRS3161261, file = 'CFTR_SRS3161261.RData')

SRS4245406 <- readMM('CFTR/Data/matrix_SRS4245406.mtx')
rownames(SRS4245406) <- readLines('CFTR/Data/geneNames_SRS4245406.txt')
colnames(SRS4245406) <- readLines('CFTR/Data/barcodes_SRS4245406.txt')
SRS4245406 <- scQC(SRS4245406)
SRS4245406 <- SRS4245406[rowMeans(SRS4245406 != 0) > 0.05,]
SRS4245406 <- SRS4245406[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(SRS4245406), ignore.case = TRUE),]
SRS4245406 <- scTenifoldKnk(SRS4245406)
save(SRS4245406, file = 'CFTR_SRS4245406.RData')
