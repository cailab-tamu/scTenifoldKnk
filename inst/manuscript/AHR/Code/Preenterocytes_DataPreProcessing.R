setwd('/data/dcosorioh/AHR/')
library(Matrix)
library(scTenifoldNet)

# WT <- readMM('Preenterocytes_WT.mtx')
# rownames(WT) <- readLines('Preenterocytes_WT_genes.txt')
# colnames(WT) <- readLines('Preenterocytes_WT_barcodes.txt')
# WT <- WT[!grepl('^MT-|^RPL|^RPS',rownames(WT), ignore.case = TRUE),]
# 
# KO <- readMM('Preenterocytes_KO.mtx')
# rownames(KO) <- readLines('Preenterocytes_KO_genes.txt')
# colnames(KO) <- readLines('Preenterocytes_KO_barcodes.txt')
# KO <- KO[!grepl('^MT-|^RPL|^RPS',rownames(KO), ignore.case = TRUE),]
# 
# O <- scTenifoldNet(WT, KO, qc_minLibSize = 0)
# save(O, file = 'dr_PreE.RData')


WT <- readMM('Preenterocytes_WT.mtx')
rownames(WT) <- readLines('Preenterocytes_WT_genes.txt')
colnames(WT) <- readLines('Preenterocytes_WT_barcodes.txt')
WT <- WT[!grepl('^MT-|^RPL|^RPS',rownames(WT), ignore.case = TRUE),]

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
  MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2)
  set.seed(1)
  DR <- scTenifoldNet::dRegulation(MA)
  outputList <- list()
  outputList$WT <- Matrix::Matrix(WT)
  outputList$KO <- Matrix::Matrix(KO)
  outputList$manifoldAlignment <- MA
  outputList$diffRegulation <- DR
  return(outputList)
}
O <- scTenifoldKnk(WT, 'Ahr')
save(O, file = 'ko_PreE.RData')
