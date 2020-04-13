source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/sccTenifoldKNK.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')

library(Matrix)

SRS3161261 <- readMM('../Data/matrix_SRS3161261.mtx')
rownames(SRS3161261) <- readLines('../Data/geneNames_SRS3161261.txt')
colnames(SRS3161261) <- readLines('../Data/barcodes_SRS3161261.txt')

SRS3161261 <- sccTenifoldKNK(SRS3161261, gKO = 'Cftr')
