library(Matrix)
library(Seurat)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')

load('Z:/Cailab/ChapkinLab_CailabReport/CailabReport_Dataset.RData')
metadata <- read.csv('Z:/Cailab/ChapkinLab_CailabReport/CailabReport_Table2.csv', stringsAsFactors = FALSE)
metadata <- metadata[metadata$CLASS == 'KO',]

ALL <- ALL@assays$RNA@counts

sapply(unique(metadata$CT), function(X){
  tempData <- ALL[,metadata$X[metadata$CT %in% X]]
  tempData <- scQC(tempData)
  if(ncol(tempData) > 500){
    tempData <- tempData[rowMeans(tempData != 0) > 0.05,]
    writeMM(tempData, paste0(X, '+KO.mtx'))
    writeLines(rownames(tempData), paste0('geneList_',X,'_KO.txt'))
    writeLines(colnames(tempData), paste0('barcodes_',X,'_KO.txt'))
  }
})
sapply(list.files(pattern = '.mtx'), function(Z){
  fileContent <- paste0("X <- '", Z,"'
         library(Matrix)
  countMatrix <- readMM(X)
  rownames(countMatrix) <- readLines(paste0('geneList_', gsub('.mtx','.txt', X)))
  colnames(countMatrix) <- readLines(paste0('barcodes_', gsub('.mtx','.txt', X)))
  
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
  O <- scTenifoldKnk(countMatrix, 'Ahr')
  save(O, file = gsub('.mtx','.RData',X))
")
  writeLines(fileContent, paste0(gsub('.mtx','.R',Z)))
  writeLines(paste0('Rscript /data/dcosorioh/AHR/',gsub('.mtx','.R',Z)),paste0(gsub('.mtx','.sh',Z)))
})
