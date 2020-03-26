#' @export scTenifoldKNK
#' @importFrom Matrix Matrix rowMeans
#' @importFrom RSpectra svds
#' @importFrom scTenifoldNet makeNetworks tensorDecomposition manifoldAlignment dRegulation
#' @importFrom rCUR CUR
#' @author Daniel Osorio <dcosorioh@tamu.edu>
#' @title scTenifoldKNK
#' @description Predict gene perturbations
#' @param countMatrix countMatrix
#' @param gKO gKO

scTenifoldKNK <- function(countMatrix, gKO = NULL, method = 'CUR', qc_mtThreshold = 0.1, qc_minLSize = 1000, nc_nNet = 10, nc_nCells = 500, nc_nComp = 3,
                          nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = 0.8){
  countMatrix <- scQC(countMatrix, mtThreshold = qc_mtThreshold, minLSize = qc_minLSize)
  countMatrix <- Matrix(vstNorm(countMatrix))
  if(ncol(countMatrix) > 500){
    countMatrix <- countMatrix[rowMeans(countMatrix != 0) >= 0.05,]
  } else {
    countMatrix[rowSums(countMatrix != 0) >= 25]
  }
  set.seed(1)
  WT <- makeNetworks(countMatrix, nNet = nc_nNet, nCells = nc_nNet, nComp = nc_nComp, scaleScores = nc_scaleScores, symmetric = nc_symmetric, q = nc_q)
  set.seed(1)
  WT <- tensorDecomposition(WT)
  WT <- as.matrix(WT$X)
  if(method == "CUR"){
    KO <- doCUR(WT, gKO)
  }
  if(method == "GGM"){

  }
  set.seed(1)
  MA <- manifoldAlignment(WT, KO)
  set.seed(1)
  DR <- dRegulation(MA, minFC = 0)
  outputList <- list()
  outputList$WT <- Matrix(WT)
  outputList$KO <- Matrix(KO)
  outputList$diffRegulation <- DR
  return(outputList)
}
