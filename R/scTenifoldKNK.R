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

scTenifoldKNK <- function(countMatrix, gKO = NULL){
  # source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/vstNorm.R')
  # source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')
  countMatrix <- scQC(countMatrix)
  countMatrix <- Matrix(vstNorm(countMatrix))
  countMatrix <- countMatrix[rowMeans(countMatrix != 0) > 0.05,]
  set.seed(1)
  WT <- makeNetworks(countMatrix)
  set.seed(1)
  WT <- tensorDecomposition(WT)
  WT <- as.matrix(WT$X)
  KO <- CUR(WT, sv = svds(WT, 5))
  C <- KO@C
  C[gKO,] <- 0
  KO <- C %*% KO@U %*% KO@R
  set.seed(1)
  MA <- manifoldAlignment(WT, KO)
  set.seed(1)
  DR <- dRegulation(MA, minFC = 0)
  outputList <- list()
  outputList$WT <- WT
  outputList$KO <- KO
  outputList$diffRegulation <- DR
  return(outputList)
}
