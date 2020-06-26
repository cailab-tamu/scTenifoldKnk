#' @export scTenifoldKNK
#' @importFrom Matrix Matrix rowMeans rowSums
#' @importFrom RSpectra svds
#' @importFrom scTenifoldNet makeNetworks tensorDecomposition manifoldAlignment dRegulation
#' @author Daniel Osorio <dcosorioh@tamu.edu>
#' @title scTenifoldKNK
#' @description Predict gene perturbations
#' @param countMatrix countMatrix
#' @param gKO gKO

sData <- read.csv('../inst/benchmark/SERGIO/simulationOutput.csv', header = FALSE)
sData <- as.matrix(sData)
rownames(sData) <- paste0('g',1:100)
colnames(sData) <- paste0('c',seq_len(ncol(sData)))

# countMatrix <- sData
scTenifoldKNK <- function(countMatrix, gKO = NULL, qc_mtThreshold = 0.1, qc_minLSize = 1000, nc_nNet = 10, nc_nCells = 500, nc_nComp = 3,
                          nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = 0.9, td_K = 3, td_maxIter = 1000,
                          td_maxError = 1e-05, ma_nDim = 2){
  #countMatrix <- scQC(countMatrix, mtThreshold = qc_mtThreshold, minLSize = qc_minLSize)
  #countMatrix <- Matrix(vstNorm(countMatrix))
  if(ncol(countMatrix) > 500){
    countMatrix <- countMatrix[rowMeans(countMatrix != 0) >= 0.05,]
  } else {
    countMatrix[rowSums(countMatrix != 0) >= 25,]
  }
  set.seed(1)
  WT <- scTenifoldNet::makeNetworks(X = countMatrix, q = nc_q, nNet = nc_nNet, nCells = nc_nCells, scaleScores = nc_scaleScores, symmetric = nc_symmetric, nComp = nc_nComp)
  set.seed(1)
  WT <- scTenifoldNet::tensorDecomposition(xList = WT, K = td_K, maxError = td_maxError, maxIter = td_maxIter)
  WT <- as.matrix(WT$X)
  WT <- (WT + t(WT))/2
  KO <- WT
  KO[gKO,] <- 0
  set.seed(1)
  MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = ma_nDim)
  set.seed(1)
  DR <- scTenifoldNet::dRegulation(MA)
  DR$FC <- (DR$distance^2)/mean(DR$distance[-seq_len(length(gKO))]^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  outputList <- list()
  outputList$WT <- Matrix(WT)
  outputList$KO <- Matrix(KO)
  outputList$MA <- MA
  outputList$diffRegulation <- DR
  return(outputList)
}

O <- scTenifoldKNK(sData, gKO = 20, nc_q = 0.9, qc_minLSize = 0)
O

source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
plotDR(O, labelGenes = 'P')
image(O$WT)
