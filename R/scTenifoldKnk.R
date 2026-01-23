#' @export scTenifoldKnk
#' @importFrom Matrix Matrix rowMeans rowSums
#' @import cli
#' @importFrom scTenifoldNet makeNetworks tensorDecomposition manifoldAlignment
#' @author Daniel Osorio <dcosorioh@tamu.edu>
#' @title scTenifoldKNK
#' @description Predict gene perturbations
#' @param countMatrix countMatrix
#' @param gKO gKO
#' @param qc A boolean value (TRUE/FALSE), if TRUE, a quality control is applied over the data.
#' @param qc_minLSize An integer value. Defines the minimum library size required for a cell to be included in the analysis.
#' @param qc_mtThreshold A decimal value between 0 and 1. Defines the maximum ratio of mitochondrial reads (mithocondrial reads / library size) present in a cell to be included in the analysis. It's computed using the symbol genes starting with 'MT-' non-case sensitive.
#' @param qc_minCells An integer value. Defines the minimum number of cells a gene should be expressed in to be included in the analysis. By default, it's set to 25.
#' @param nc_nNet An integer value. The number of networks based on principal components regression to generate.
#' @param nc_nCells An integer value. The number of cells to subsample each time to generate a network.
#' @param nc_nComp An integer value. The number of principal components in PCA to generate the networks. Should be greater than 2 and lower than the total number of genes.
#' @param nc_symmetric A boolean value (TRUE/FALSE), if TRUE, the weights matrix returned will be symmetric.
#' @param nc_scaleScores A boolean value (TRUE/FALSE), if TRUE, the weights will be normalized such that the maximum absolute value is 1.
#' @param nc_lambda A continuous value between 0 and 1. Defines the multiplicative value (1-lambda) to be applied over the weaker edge connecting two genes to maximize the adjacency matrix directionality.
#' @param nc_q A decimal value between 0 and 1. Defines the cut-off threshold of top q\% relationships to be returned.
#' @param td_K An integer value. Defines the number of rank-one tensors used to approximate the data using CANDECOMP/PARAFAC (CP) Tensor Decomposition.
#' @param td_maxIter An integer value. Defines the maximum number of iterations if error stay above \code{td_maxError}.
#' @param td_maxError A decimal value between 0 and 1. Defines the relative Frobenius norm error tolerance.
#' @param td_nDecimal An integer value indicating the number of decimal places to be used.
#' @param ma_nDim An integer value. Defines the number of dimensions of the low-dimensional feature space to be returned from the non-linear manifold alignment.
#' @param nCores An integer value. Defines the number of cores to be used.
#' @examples
#' # Loading single-cell data
#' scRNAseq <- system.file("single-cell/example.csv", package = "scTenifoldKnk")
#' scRNAseq <- read.csv(scRNAseq, row.names = 1)
#'
#' # Running scTenifoldKnk
#' scTenifoldKnk(countMatrix = scRNAseq, gKO = "G100", qc_minLSize = 0)
scTenifoldKnk <- function(countMatrix, qc = TRUE, gKO = NULL, qc_mtThreshold = 0.1, qc_minLSize = 1000, qc_minCells = 25, nc_lambda = 0, nc_nNet = 10, nc_nCells = 500, nc_nComp = 3,
                          nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = 0.9, td_K = 3, td_maxIter = 1000,
                          td_maxError = 1e-05, td_nDecimal = 3, ma_nDim = 2, nCores = parallel::detectCores()) {
  
  # Start a CLI process to report progress to the user
  cli::cli_h1("scTenifoldKnk")
  cli::cli_alert_info("Simulating {gKO} gene knockout")

  # Check that the requested gene to knock out is present in the input matrix
  if (!gKO %in% rownames(countMatrix)) {
    cli::cli_alert_danger("{gKO} is not present in the count matrix used as input")
    cli::cli_abort("{gKO} is not present in the count matrix used as input")
  }

  # Optional quality control: filter cells and genes
  if (isTRUE(qc)) {
    # Remove low-quality cells (using mitochondrial ratio and library size)
    countMatrix <- scQC(countMatrix, mtThreshold = qc_mtThreshold, minLSize = qc_minLSize)
    # Keep genes expressed in at least `qc_minCells` cells (row subset)
    countMatrix[rowSums(countMatrix != 0) >= qc_minCells, ]
    # Report QC results
    cli::cli_alert_success("Count matrix quality control applied: retained {nrow(countMatrix)} genes and {ncol(countMatrix)} cells")
  }

  # Re-check presence of the KO gene after filtering
  if (!gKO %in% rownames(countMatrix)) {
    cli::cli_alert_danger("{gKO} is not present in the count matrix after count matrix quality control")
    cli::cli_abort("{gKO} is not present in the count matrix after count matrix quality control")
  }

  # Build an ensemble of gene regulatory networks (subsample cells, use PCR)
  WT <- scTenifoldNet::makeNetworks(X = countMatrix, q = nc_q, nNet = nc_nNet, nCells = nc_nCells, scaleScores = nc_scaleScores, symmetric = nc_symmetric, nComp = nc_nComp, nCores = nCores)
  cli::cli_alert_success("Network construction complete (nNet = {nc_nNet}, nCells per net = {nc_nCells})")

  # Tensor decomposition (CP) to denoise / approximate the ensemble of networks
  WT <- scTenifoldNet::tensorDecomposition(xList = WT, K = td_K, maxError = td_maxError, maxIter = td_maxIter, nDecimal = td_nDecimal)
  cli::cli_alert_success("Tensor decomposition completed (K = {td_K})")

  # Extract reconstructed network and enforce directionality
  WT <- WT$X
  WT <- strictDirection(WT, lambda = nc_lambda)
  WT <- as.matrix(WT)
  # Remove self-loops
  diag(WT) <- 0
  # Transpose to have genes as rows for downstream steps
  WT <- t(WT)
  cli::cli_alert_success("Prepared WT adjacency matrix for KO simulation")

  # Simulate knockout by zeroing outgoing edges from the KO gene
  KO <- WT
  KO[gKO, ] <- 0
  cli::cli_alert_success("Simulated knockout for {gKO}")

  # Align WT and KO networks into a shared low-dimensional manifold space
  MA <- manifoldAlignment(WT, KO, d = ma_nDim, nCores = nCores)
  cli::cli_alert_success("Manifold alignment completed (d = {ma_nDim})")

  # Compute differential regulation for the KO gene from the aligned manifolds
  DR <- dRegulation(MA, gKO)
  cli::cli_alert_success("Differential regulation computed for {gKO}")

  # Prepare and return results
  outputList <- list()
  outputList$tensorNetworks$WT <- Matrix(WT)
  outputList$tensorNetworks$KO <- Matrix(KO)
  outputList$manifoldAlignment <- MA
  outputList$diffRegulation <- DR
  # Finish CLI process and return results
  cli::cli_alert_success("Finished scTenifoldKnk for {gKO}")
  return(outputList)
}
