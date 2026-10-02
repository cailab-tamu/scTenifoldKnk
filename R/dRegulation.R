#' @export dRegulation
#' @importFrom stats dist pchisq pnorm p.adjust qqnorm
#' @importFrom MASS boxcox
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_warning
#' @title Evaluates gene differential regulation based on manifold alignment
#'   distances.
#' @description Using the output of the non-linear manifold alignment, this
#'   function computes the Euclidean distance between the coordinates for the
#'   same gene in both conditions. Calculated distances are then transformed
#'   using Box-Cox power transformation, and standardized to ensure normality.
#'   P-values are assigned following the chi-square distribution over the
#'   fold-change of the squared distance computed with respect to the
#'   expectation, the mean squared distance of the genes that were not knocked
#'   out (\code{gKO}), or, when \code{empiricalNull = TRUE}, using Efron's
#'   empirical null estimated from the Z-scores with \code{locfdr}. Genes whose
#'   distance is at the level of floating-point noise (at most
#'   \code{sqrt(.Machine$double.eps)} times the largest absolute coordinate) did
#'   not move between conditions and get a p-value of 1; if this applies to
#'   every gene, for example when the knockout leaves the network unchanged, a
#'   warning is raised.
#' @param manifoldOutput A matrix. The output of the non-linear manifold
#'   alignment, a labeled matrix with two times the number of shared genes as
#'   rows (X_ genes followed by Y_ genes in the same order) and \code{d} number
#'   of columns.
#' @param gKO A character vector with the knocked-out genes, or \code{NULL}.
#'   These genes are perturbed by construction, so they are left out of the
#'   expectation used to compute the fold-changes; otherwise their large
#'   distances inflate it and hide the other genes. They are still reported in
#'   the output. Default: NULL, the expectation is computed from all genes.
#' @param empiricalNull A boolean value (TRUE/FALSE). If TRUE, p-values are
#'   assigned using Efron's empirical null: the null distribution of the
#'   Z-scores is estimated from the bulk of the data with \code{locfdr::locfdr}
#'   instead of assuming the theoretical chi-square null. Requires the
#'   \code{locfdr} package. Default: FALSE.
#' @param direction Optional named numeric vector of direction scores (names:
#'   genes), for example a row of \code{\link{knockoutDirection}}. If given,
#'   the columns \code{direction} and \code{directionScore} are appended.
#'   Default: NULL.
#' @return A data frame with 6 columns (8 when \code{direction} is given) as follows: \itemize{
#' \item \code{gene} A character vector with the gene id identified from the
#'   \code{manifoldAlignment} output.
#' \item \code{distance} A numeric vector of the Euclidean distance computed
#'   between the coordinates of the same gene in both conditions.
#' \item \code{Z} A numeric vector of the Z-scores computed after Box-Cox power
#'   transformation.
#' \item \code{FC} A numeric vector of the FC computed with respect to the
#'   expectation, the mean squared distance of the genes not in \code{gKO}.
#' \item \code{p.value} A numeric vector of the p-values associated to the
#'   fold-changes, probabilities are assigned as \eqn{P[X > x]} using the
#'   Chi-square distribution with one degree of freedom.
#' \item \code{p.adj} A numeric vector of adjusted p-values using Benjamini &
#'   Hochberg (1995) FDR correction.
#' \item \code{direction} Only when \code{direction} is given: \code{"up"} or
#'   \code{"down"}, the predicted direction of the change of the gene; the
#'   knocked-out genes are \code{"down"} by construction.
#' \item \code{directionScore} Only when \code{direction} is given: the
#'   signed direction score (\code{NA} for the knocked-out genes).
#' }
#' @references \itemize{
#' \item Benjamini, Y., and Yekutieli, D. (2001). The control of the false
#'   discovery rate in multiple testing under dependency. Annals of Statistics,
#'   29, 1165-1188. doi: 10.1214/aos/1013699998.
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
#' # Performing single-cell quality control
#' qcOutput <- scQC(
#'   X = X,
#'   minLibSize = 30,
#'   removeOutlierCells = TRUE,
#'   minPCT = 0.05,
#'   maxMTratio = 0.1
#' )
#'
#' # Computing 3 gene regulatory networks from subsamples of 500 cells
#' xNetworks <- scTenifoldNet::makeNetworks(
#'   X = qcOutput,
#'   nNet = 3,
#'   nCells = 500,
#'   nComp = 3,
#'   scaleScores = TRUE,
#'   symmetric = FALSE,
#'   q = 0.95
#' )
#'
#' # Computing a K = 3 CANDECOMP/PARAFAC (CP) Tensor Decomposition
#' tdOutput <- scTenifoldNet::tensorDecomposition(xNetworks, K = 3, maxError = 1e5, maxIter = 1e3)
#'
#' \dontrun{
#' # Computing the manifold alignment
#' maOutput <- scTenifoldNet::manifoldAlignment(tdOutput$X, tdOutput$X)
#'
#' # Evaluating differential regulation
#' drOutput <- dRegulation(maOutput)
#' head(drOutput)
#'
#' # Plotting — genes with FDR < 0.05 colored in red
#' geneColor <- ifelse(drOutput$p.adj < 0.05, 'red', 'black')
#' qqnorm(drOutput$Z, main = 'Standardized Distance', pch = 16, col = geneColor)
#' qqline(drOutput$Z)
#' }

dRegulation <- function(manifoldOutput, gKO = NULL, empiricalNull = FALSE, direction = NULL) {

  geneList <- rownames(manifoldOutput)
  geneList <- geneList[grepl('^X_', geneList)]
  geneList <- gsub('^X_', '', geneList)
  nGenes <- length(geneList)

  eGenes <- nrow(manifoldOutput) / 2

  eGeneList <- rownames(manifoldOutput)
  eGeneList <- eGeneList[grepl('^Y_', eGeneList)]
  eGeneList <- gsub('^Y_', '', eGeneList)

  if (nGenes != eGenes) {
    stop('Number of identified and expected genes are not the same')
  }
  if (!all(eGeneList == geneList)) {
    stop('Genes are not ordered as expected. ',
         'X_ genes should be followed by Y_ genes in the same order')
  }

  if (!is.null(gKO) && (!is.character(gKO) || anyNA(gKO))) {
    stop("'gKO' must be a character vector of gene symbols")
  }
  missingGenes <- setdiff(gKO, geneList)
  if (length(missingGenes) > 0) {
    stop("The following genes are not present in the manifold alignment: ",
         paste(missingGenes, collapse = ", "))
  }
  isKO <- geneList %in% gKO
  if (all(isKO)) {
    stop("At least one gene that was not knocked out is required to compute the expectation")
  }

  cli::cli_alert_info("Computing distances for {nGenes} genes")

  dMetric <- sapply(seq_len(nGenes), function(G) {
    X <- manifoldOutput[G, ]
    Y <- manifoldOutput[(G + nGenes), ]
    as.numeric(dist(rbind(X, Y)))
  })

  # Distances at the level of floating-point noise mean the gene did not move
  # between conditions; ranking them against each other would flag noise
  noiseLevel <- sqrt(.Machine$double.eps) * max(abs(manifoldOutput))
  .drStatistics(dMetric, geneList, isKO, empiricalNull, noiseLevel, direction)
}

# Box-Cox / Z-score / chi-square statistics of the per-gene distances, shared by
# the manifold alignment and the heat manifold alignment routes
.drStatistics <- function(dMetric, geneList, isKO, empiricalNull = FALSE, noiseLevel = 0,
                          direction = NULL) {
  nGenes <- length(geneList)
  # Box-Cox transformation
  lambdaValues <- seq(-2, 2, length.out = 1000)
  lambdaValues <- lambdaValues[lambdaValues != 0]
  BC <- try(MASS::boxcox(dMetric ~ 1, plot = FALSE, lambda = lambdaValues),
            silent = TRUE)
  if (inherits(BC, 'try-error')) {
    nD <- dMetric
  } else {
    BC <- BC$x[which.max(BC$y)]
    if (BC < 0) {
      nD <- 1 / (dMetric ^ BC)
    } else {
      nD <- dMetric ^ BC
    }
  }

  Z <- scale(nD)
  # The knocked-out genes are perturbed by construction; their large
  # distances would inflate the expectation and hide the other genes
  E <- mean(dMetric[!isKO]^2)
  FC <- dMetric^2 / E

  if (isTRUE(empiricalNull)) {
    if (!requireNamespace('locfdr', quietly = TRUE)) {
      stop("Package 'locfdr' is required when 'empiricalNull = TRUE'. ",
           "Install it with install.packages('locfdr').")
    }
    zScores <- as.numeric(Z)
    eNull <- try(locfdr::locfdr(zScores, plot = 0), silent = TRUE)
    if (inherits(eNull, 'try-error')) {
      cli::cli_alert_warning(
        "Empirical null estimation failed; falling back to the theoretical null"
      )
      pValues <- pchisq(q = FC, df = 1, lower.tail = FALSE)
    } else {
      delta0 <- eNull$fp0['mlest', 'delta']
      sigma0 <- eNull$fp0['mlest', 'sigma']
      cli::cli_alert_info(
        "Efron empirical null: mean = {round(delta0, 3)}, sd = {round(sigma0, 3)}"
      )
      pValues <- pnorm((zScores - delta0) / sigma0, lower.tail = FALSE)
    }
  } else {
    pValues <- pchisq(q = FC, df = 1, lower.tail = FALSE)
  }

  isNoise <- dMetric <= noiseLevel
  pValues[isNoise] <- 1
  if (all(isNoise)) {
    warning("No gene differs between the two conditions beyond numerical noise; ",
            "all p-values were set to 1")
  }
  pAdjusted <- p.adjust(pValues, method = 'fdr')

  dOut <- data.frame(
    gene = geneList,
    distance = dMetric,
    Z = Z,
    FC = FC,
    p.value = pValues,
    p.adj = pAdjusted
  )
  if (!is.null(direction)) {
    if (!is.numeric(direction) || is.null(names(direction))) {
      stop("'direction' must be a named numeric vector of direction scores")
    }
    score <- unname(direction[geneList])
    score[isKO] <- NA
    dOut$direction <- ifelse(isKO | score < 0, "down", ifelse(score > 0, "up", NA_character_))
    dOut$directionScore <- score
  }
  dOut <- dOut[order(dOut$p.value), ]
  dOut <- as.data.frame.array(dOut)

  nSig <- sum(dOut$p.adj < 0.05)
  cli::cli_alert_success(
    "Differential regulation complete: {nSig}/{nGenes} significant genes (FDR < 0.05)"
  )
  return(dOut)
}
