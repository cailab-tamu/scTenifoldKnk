#' @export scTenifoldKnk
#' @importFrom methods as
#' @importFrom cli cli_h1 cli_alert_info cli_alert_success
#' @importFrom Matrix Matrix
#' @importFrom scTenifoldNet makeNetworks tensorDecomposition manifoldAlignment cpmNormalization
#' @author Daniel Osorio <dcosorioh@gmail.com>
#' @title scTenifoldKNK
#' @description Predict gene perturbations using in-silico knockout experiments
#'   from single-cell gene regulatory networks.
#' @param countMatrix Raw counts matrix with cells as columns and genes (symbols) as rows, as a \code{matrix} or a sparse \code{dgCMatrix}. A data.frame is not accepted; convert it with \code{as.matrix()} first.
#' @param gKO Character. In knockout mode (\code{transcriptomeWide = FALSE}), the gene symbol of the gene to knock out, or a character vector of several genes to knock out together in a single simulated experiment (e.g. \code{c("Hnf4a", "Hnf4g")}). In transcriptome-wide mode (\code{transcriptomeWide = TRUE}), an optional character vector defining the subset of genes to perturb, each one knocked out separately; if \code{NULL}, every gene in the WT network is perturbed.
#' @param transcriptomeWide A boolean value (TRUE/FALSE). If TRUE, the WT network is built once and each target gene is knocked out in turn, returning the manifold-alignment distances for every perturbation. Default: FALSE.
#' @param qc A boolean value (TRUE/FALSE), if TRUE, a quality control is applied over the data.
#' @param qc_minLibSize An integer value. Defines the minimum library size required for a cell to be included in the analysis.
#' @param qc_removeOutlierCells A boolean value (TRUE/FALSE), if TRUE, cells with library size identified as outliers are removed. For further details see: \code{?boxplot.stats}
#' @param qc_minPCT A decimal value between 0 and 1. Defines the minimum fraction of cells where the gene needs to be expressed to be included in the analysis.
#' @param qc_maxMTratio A decimal value between 0 and 1. Defines the maximum ratio of mitochondrial reads (mitochondrial reads / library size) present in a cell to be included in the analysis. It's computed using the symbol genes starting with 'MT-' non-case sensitive.
#' @param nc_nNet An integer value. The number of networks based on principal components regression to generate.
#' @param nc_nCells An integer value. The number of cells to subsample each time to generate a network.
#' @param nc_nComp An integer value. The number of principal components in PCA to generate the networks. Should be greater than 2 and lower than the total number of genes.
#' @param nc_symmetric A boolean value (TRUE/FALSE), if TRUE, the weights matrix returned will be symmetric.
#' @param nc_scaleScores A boolean value (TRUE/FALSE), if TRUE, the weights will be normalized such that the maximum absolute value is 1.
#' @param nc_lambda A continuous value between 0 and 1. Defines the multiplicative value (1-lambda) to be applied over the weaker edge connecting two genes to maximize the adjacency matrix directionality.
#' @param nc_q A decimal value between 0 and 1. Defines the cut-off threshold of top q\% relationships to be returned.
#' @param nc_priorNetwork A data.frame containing a prior gene regulatory network. The data.frame must have two columns: `regulators` and `targets`. Default: NULL.
#' @param td_K An integer value. Defines the number of rank-one tensors used to approximate the data using CANDECOMP/PARAFAC (CP) Tensor Decomposition.
#' @param td_maxIter An integer value. Defines the maximum number of iterations if error stay above \code{td_maxError}.
#' @param td_maxError A decimal value between 0 and 1. Defines the relative Frobenius norm error tolerance.
#' @param td_nDecimal An integer value indicating the number of decimal places to be used.
#' @param ma_nDim An integer value. Defines the number of dimensions of the low-dimensional feature space to be returned from the non-linear manifold alignment.
#' @param dr_empiricalNull A boolean value (TRUE/FALSE). If TRUE, the differential regulation p-values are assigned using Efron's empirical null (estimated with \code{locfdr}) instead of the theoretical chi-square null. Requires the \code{locfdr} package. Default: FALSE.
#' @param nCores An integer value. Defines the number of cores to be used.
#' @param seed An integer value. The RNG is set to this seed before each random stage (network construction, tensor decomposition and manifold alignment), so results are reproducible and independent of the caller's RNG state; the caller's RNG state is restored on exit. Use different values to assess run-to-run variability. If \code{NULL}, the RNG is never reseeded and the caller's RNG state (e.g. a previous \code{set.seed()}) drives all random stages. Default: 1.
#' @return In single knockout mode (\code{transcriptomeWide = FALSE}), a list with 3 slots as follows:
#' \itemize{
#' \item{tensorNetworks:} The WT and KO weight-averaged denoised gene regulatory networks.
#' \item{manifoldAlignment:} The generated low-dimensional features result of the non-linear manifold alignment.
#' \item{diffRegulation:} The results of the differential regulation analysis.
#' }
#' In transcriptome-wide mode (\code{transcriptomeWide = TRUE}), a list with 2 slots as follows:
#' \itemize{
#' \item{tensorNetworks:} A list with the WT weight-averaged denoised gene regulatory network.
#' \item{perturbationDistances:} A numeric matrix of manifold-alignment distances with the perturbed genes as rows and all genes in the WT network as columns.
#' }
#' @examples
#' library(scTenifoldKnk)
#'
#' # Simulating a dataset following a negative binomial distribution with high sparsity (~67%)
#' nCells = 2000
#' nGenes = 100
#' set.seed(1)
#' X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
#' X <- round(X)
#' X <- matrix(X, ncol = nCells)
#' rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
#'
#' \dontrun{
#' # Running scTenifoldKnk — simulating knockout of gene ng10
#' output <- scTenifoldKnk(
#'   countMatrix = X,
#'   gKO = "ng10",
#'   nc_nNet = 10,
#'   nc_nCells = 500,
#'   td_K = 3,
#'   qc_minLibSize = 30
#' )
#'
#' # Structure of the output
#' str(output)
#'
#' # Accessing the WT and KO gene regulatory networks
#' dim(output$tensorNetworks$WT)
#' dim(output$tensorNetworks$KO)
#'
#' # Accessing the manifold alignment result
#' head(output$manifoldAlignment)
#'
#' # Differential regulation results — top perturbed genes
#' head(output$diffRegulation, n = 10)
#'
#' # Plotting the KO-centered subnetwork
#' plotKO(output, gKO = "ng10")
#'
#' # Multi-gene knockout — ng10 and ng20 knocked out together
#' dkoOutput <- scTenifoldKnk(
#'   countMatrix = X,
#'   gKO = c("ng10", "ng20"),
#'   nc_nNet = 10,
#'   nc_nCells = 500,
#'   td_K = 3,
#'   qc_minLibSize = 30
#' )
#' head(dkoOutput$diffRegulation, n = 10)
#' plotKO(dkoOutput, gKO = c("ng10", "ng20"))
#'
#' # Transcriptome-wide perturbation — knock out every gene in the WT network
#' twOutput <- scTenifoldKnk(
#'   countMatrix = X,
#'   transcriptomeWide = TRUE,
#'   nc_nNet = 10,
#'   nc_nCells = 500,
#'   td_K = 3,
#'   qc_minLibSize = 30
#' )
#'
#' # Matrix of distances: perturbed genes (rows) by all genes (columns)
#' dim(twOutput$perturbationDistances)
#' twOutput$perturbationDistances[1:5, 1:5]
#' }
scTenifoldKnk <- function(countMatrix, gKO = NULL, transcriptomeWide = FALSE,
                          qc = TRUE,
                          qc_minLibSize = 1000, qc_removeOutlierCells = TRUE,
                          qc_minPCT = 0.05, qc_maxMTratio = 0.1,
                          nc_lambda = 0, nc_nNet = 10, nc_nCells = 500,
                          nc_nComp = 3, nc_scaleScores = TRUE,
                          nc_symmetric = FALSE, nc_q = 0.9,
                          nc_priorNetwork = NULL, td_K = 3,
                          td_maxIter = 1000, td_maxError = 1e-05,
                          td_nDecimal = 3, ma_nDim = 2,
                          dr_empiricalNull = FALSE,
                          nCores = parallel::detectCores(),
                          seed = 1) {

  cli::cli_h1("scTenifoldKnk Pipeline")

  if (!is.null(seed)) {
    oldSeed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(restoreSeed(oldSeed), add = TRUE)
  }

  if (is.data.frame(countMatrix)) {
    stop("'countMatrix' must be a matrix or a sparse dgCMatrix, not a data.frame. ",
         "Convert it first with countMatrix <- as.matrix(countMatrix)")
  }

  # gKO is optional in transcriptome-wide mode (every gene is perturbed) and
  # required otherwise; one or more genes are knocked out together
  if (!is.null(gKO) && (!is.character(gKO) || length(gKO) == 0 || anyNA(gKO))) {
    stop("'gKO' must be a character vector of gene symbols")
  }
  if (!isTRUE(transcriptomeWide) && is.null(gKO)) {
    stop("At least one gene symbol must be provided in 'gKO' to perform the knockout")
  }
  gKO <- unique(gKO)

  # Check that the requested genes are present in the input matrix
  missingGenes <- gKO[!gKO %in% rownames(countMatrix)]
  if (length(missingGenes) > 0) {
    stop("The following genes are not present in the count matrix used as input: ",
         paste(missingGenes, collapse = ", "))
  }

  # Step 1: Quality Control
  if (isTRUE(qc)) {
    cli::cli_alert_info("Step 1/7: Quality control")
    countMatrix <- scQC(countMatrix, minLibSize = qc_minLibSize,
                        removeOutlierCells = qc_removeOutlierCells,
                        minPCT = qc_minPCT, maxMTratio = qc_maxMTratio)
  }

  # Re-check presence of the KO gene(s) after filtering
  missingGenes <- gKO[!gKO %in% rownames(countMatrix)]
  if (length(missingGenes) > 0) {
    stop("The following genes are not present in the count matrix after quality control: ",
         paste(missingGenes, collapse = ", "))
  }

  # Step 2: CPM Normalization
  cli::cli_alert_info("Step 2/7: CPM normalization")
  countMatrix <- cpmNormalization(countMatrix)

  # Step 3: Network construction
  cli::cli_alert_info("Step 3/7: Building gene regulatory networks")
  if (!is.null(seed)) set.seed(seed)
  WT <- makeNetworks(X = countMatrix, q = nc_q,
                     priorNetwork = nc_priorNetwork, nNet = nc_nNet,
                     nCells = nc_nCells, scaleScores = nc_scaleScores,
                     symmetric = nc_symmetric, nComp = nc_nComp,
                     nCores = nCores)

  # Step 4: Tensor decomposition
  cli::cli_alert_info("Step 4/7: Tensor decomposition")
  if (!is.null(seed)) set.seed(seed)
  WT <- tensorDecomposition(xList = WT, K = td_K, maxError = td_maxError,
                            maxIter = td_maxIter, nDecimal = td_nDecimal,
                            seed = seed)

  # Extract reconstructed network, enforce directionality
  WT <- WT$X
  WT <- strictDirection(WT, lambda = nc_lambda)
  WT <- as.matrix(WT)
  diag(WT) <- 0
  WT <- t(WT)

  if (isTRUE(transcriptomeWide)) {
    # Transcriptome-wide mode: perturb each target gene in the WT network
    geneList <- rownames(WT)
    targetGenes <- if (is.null(gKO)) geneList else gKO

    # Some subset genes may be dropped during network construction
    missingGenes <- targetGenes[!targetGenes %in% geneList]
    if (length(missingGenes) > 0) {
      stop("The following genes are not present in the WT network: ",
           paste(missingGenes, collapse = ", "))
    }

    cli::cli_alert_info(
      "Step 5/5: Perturbing {length(targetGenes)} gene{?s} transcriptome-wide"
    )

    perturbationDistances <- matrix(
      NA_real_, nrow = length(targetGenes), ncol = length(geneList),
      dimnames = list(targetGenes, geneList)
    )

    cli::cli_progress_bar("Perturbing genes", total = length(targetGenes))
    for (g in targetGenes) {
      KO <- WT
      KO[g, ] <- 0
      if (!is.null(seed)) set.seed(seed)
      MA <- manifoldAlignment(WT, KO, d = ma_nDim, nCores = nCores)
      DR <- dRegulation(MA, empiricalNull = dr_empiricalNull)
      perturbationDistances[g, DR$gene] <- DR$distance
      cli::cli_progress_update()
    }
    cli::cli_progress_done()

    outputList <- list()
    outputList$tensorNetworks$WT <- Matrix(WT)
    outputList$perturbationDistances <- perturbationDistances

    cli::cli_alert_success("scTenifoldKnk pipeline complete")
    return(outputList)
  }

  # Step 5: Simulate knockout by zeroing outgoing edges from the KO gene(s)
  cli::cli_alert_info("Step 5/7: Simulating {gKO} knockout")
  KO <- WT
  KO[gKO, ] <- 0

  # A gene without outgoing edges leaves the network unchanged, so any hits are noise
  noEdges <- gKO[rowSums(WT[gKO, , drop = FALSE] != 0) == 0]
  if (length(noEdges) == length(gKO)) {
    warning(paste(gKO, collapse = ", "),
            if (length(gKO) == 1) " has" else " have",
            " no outgoing edges in the WT network; the knockout does not ",
            "change the network and the differential regulation results reflect ",
            "only numerical noise. Consider 'dr_empiricalNull = TRUE'.")
  } else if (length(noEdges) > 0) {
    warning("The following genes have no outgoing edges in the WT network, so ",
            "knocking them out has no effect: ", paste(noEdges, collapse = ", "))
  }

  # Step 6: Manifold alignment
  cli::cli_alert_info("Step 6/7: Manifold alignment")
  if (!is.null(seed)) set.seed(seed)
  MA <- manifoldAlignment(WT, KO, d = ma_nDim, nCores = nCores)

  # Step 7: Differential regulation analysis
  cli::cli_alert_info("Step 7/7: Differential regulation analysis")
  DR <- dRegulation(MA, empiricalNull = dr_empiricalNull)

  outputList <- list()
  outputList$tensorNetworks$WT <- Matrix(WT)
  outputList$tensorNetworks$KO <- Matrix(KO)
  outputList$manifoldAlignment <- MA
  outputList$diffRegulation <- DR

  cli::cli_alert_success("scTenifoldKnk pipeline complete")
  return(outputList)
}
