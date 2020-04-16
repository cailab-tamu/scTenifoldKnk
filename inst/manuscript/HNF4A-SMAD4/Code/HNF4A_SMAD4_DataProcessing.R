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

GSM3477499 <- Seurat::Read10X_h5('HNF4A-HNF4G/Data/GSM3477499_WT_ScRNAseq_filtered_gene_bc_matrices.h5')
GSM3477499 <- scQC(GSM3477499, mtThreshold = 1)
GSM3477499 <- GSM3477499[rowMeans(GSM3477499 != 0) > 0.05,]
GSM3477499 <- GSM3477499[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(GSM3477499), ignore.case = TRUE),]
GSM3477499 <- scTenifoldKnk(GSM3477499, gKO = c('Hnf4a','Smad4'))
save(GSM3477499, file = 'HNF4A_SMAD4_GSM4116571.RData')

