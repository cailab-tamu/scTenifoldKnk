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

WT <- readMM('NKX2-1/Data/GSM3716703_Nkx2-1_control_scRNAseq_matrix.mtx.gz')
rownames(WT) <- read.csv('NKX2-1/Data/GSM3716703_Nkx2-1_control_scRNAseq_genes.tsv.gz', sep = '\t', header = FALSE, stringsAsFactors = FALSE)[,2]
colnames(WT) <- readLines('NKX2-1/Data/GSM3716703_Nkx2-1_control_scRNAseq_barcodes.tsv.gz')
WT <- scQC(WT)
WT <- WT[rowMeans(WT != 0) > 0.05,]
WT <- WT[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(WT), ignore.case = TRUE),]
WT <- scTenifoldKnk(WT, gKO = 'Nkx2-1')
save(WT, file = 'NKX21_GSM3716703.RData')

