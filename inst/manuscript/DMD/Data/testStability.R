args = commandArgs(trailingOnly=TRUE)

setwd('/data/dcosorioh/knkStability')
library(scTenifoldKnk)
library(Matrix)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')

X <- readMM('GSM4116571_Qc_matrix.mtx.gz')
colnames(X) <- readLines('GSM4116571_Qc_barcodes.tsv.gz')
rownames(X) <- read.csv('GSM4116571_Qc_features.tsv.gz', header = FALSE, sep = '\t')[,2]
X <- scQC(X)

tempX <- X[,sample(seq_len(ncol(X)),1000)]
tempX <- tempX[rowMeans(X != 0) > 0.05,]
oX <- scTenifoldKnk(countMatrix = tempX, gKO = 'Dmd')
write.csv(oX$diffRegulation, file = paste0('dmdStability_',args[1], '.csv'), row.names = FALSE)
