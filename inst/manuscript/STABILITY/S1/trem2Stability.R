args = commandArgs(trailingOnly=TRUE)

setwd('/data/dcosorioh/knkStability')
library(scTenifoldKnk)
library(Matrix)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')

X <- readMM('matrix_WT.mtx')
rownames(X) <- readLines('geneNames_WT.txt')
colnames(X) <- readLines('barcodes_WT.txt')

tempX <- X[,sample(seq_len(ncol(X)),500)]
tempX <- tempX[rowMeans(X != 0) > 0.05,]
oX <- scTenifoldKnk(countMatrix = tempX, gKO = 'Trem2')
write.csv(oX$diffRegulation, file = paste0('trem2Stability_',args[1], '.csv'), row.names = FALSE)
