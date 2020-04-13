# source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/sccTenifoldKNK.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')

library(Matrix)

SRS3161261 <- readMM('../Data/matrix_SRS3161261.mtx')
rownames(SRS3161261) <- readLines('../Data/geneNames_SRS3161261.txt')
colnames(SRS3161261) <- readLines('../Data/barcodes_SRS3161261.txt')

SRS3161261 <- sccTenifoldKNK(SRS3161261, gKO = 'Cftr')
# SRS3161261$KO <- SRS3161261$WT
# SRS3161261$KO['Cftr',] <- 0
# SRS3161261$manifoldAlignment <- scTenifoldNet::manifoldAlignment(SRS3161261$WT, SRS3161261$KO)
SRS3161261$diffRegulation <- scTenifoldNet::dRegulation(SRS3161261$manifoldAlignment[,1:2])
DR <- SRS3161261$diffRegulation
plotDR(SRS3161261)
plotKO(SRS3161261, 'Cftr')
