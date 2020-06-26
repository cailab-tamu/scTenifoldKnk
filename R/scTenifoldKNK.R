#' @export scTenifoldKNK
#' @importFrom Matrix Matrix rowMeans rowSums
#' @importFrom RSpectra svds
#' @importFrom scTenifoldNet makeNetworks tensorDecomposition manifoldAlignment dRegulation
#' @author Daniel Osorio <dcosorioh@tamu.edu>
#' @title scTenifoldKNK
#' @description Predict gene perturbations
#' @param countMatrix countMatrix
#' @param gKO gKO
#' @param qc_minLSize An integer value. Defines the minimum library size required for a cell to be included in the analysis.
#' @param qc_mtThreshold A decimal value between 0 and 1. Defines the maximum ratio of mitochondrial reads (mithocondrial reads / library size) present in a cell to be included in the analysis. It's computed using the symbol genes starting with 'MT-' non-case sensitive.
#' @param nc_nNet An integer value. The number of networks based on principal components regression to generate.
#' @param nc_nCells An integer value. The number of cells to subsample each time to generate a network.
#' @param nc_nComp An integer value. The number of principal components in PCA to generate the networks. Should be greater than 2 and lower than the total number of genes.
#' @param nc_symmetric A boolean value (TRUE/FALSE), if TRUE, the weights matrix returned will be symmetric.
#' @param nc_scaleScores A boolean value (TRUE/FALSE), if TRUE, the weights will be normalized such that the maximum absolute value is 1.
#' @param nc_q A decimal value between 0 and 1. Defines the cut-off threshold of top q\% relationships to be returned.
#' @param td_K An integer value. Defines the number of rank-one tensors used to approximate the data using CANDECOMP/PARAFAC (CP) Tensor Decomposition. 
#' @param td_maxIter An integer value. Defines the maximum number of iterations if error stay above \code{td_maxError}.
#' @param td_maxError A decimal value between 0 and 1. Defines the relative Frobenius norm error tolerance.
#' @param ma_nDim An integer value. Defines the number of dimensions of the low-dimensional feature space to be returned from the non-linear manifold alignment.

scTenifoldKnk <- function(countMatrix, gKO = NULL, qc_mtThreshold = 0.1, qc_minLSize = 1000, nc_nNet = 10, nc_nCells = 500, nc_nComp = 3,
                          nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = 0.9, td_K = 3, td_maxIter = 1000,
                          td_maxError = 1e-05, ma_nDim = 2){
  countMatrix <- scQC(countMatrix, mtThreshold = qc_mtThreshold, minLSize = qc_minLSize)
  countMatrix <- Matrix(vstNorm(countMatrix))
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
  KO <- WT
  KO[gKO,] <- 0
  set.seed(1)
  MA <- manifoldAlignment(WT, KO, d = ma_nDim)
  set.seed(1)
  DR <- dRegulation(MA, gKO)
  outputList <- list()
  outputList$WT <- Matrix(WT)
  outputList$KO <- Matrix(KO)
  outputList$MA <- MA
  outputList$diffRegulation <- DR
  return(outputList)
}