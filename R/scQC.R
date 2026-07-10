#' @export scQC
#' @title Performs single-cell data quality control
#' @importFrom grDevices boxplot.stats
#' @importFrom Matrix colSums rowMeans
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_warning
#' @description This function performs quality control filters over the provided
#'   input matrix. It checks for minimum cell library size, mitochondrial ratio,
#'   outlier cells, and the fraction of cells where a gene is expressed.
#' @param X Raw counts matrix with cells as columns and genes (symbols) as rows.
#' @param minLibSize An integer value. Defines the minimum library size required
#'   for a cell to be included in the analysis.
#' @param removeOutlierCells A boolean value (\code{TRUE/FALSE}), if
#'   \code{TRUE}, the identified cells with library size greater than 1.58
#'   IQR/sqrt(n) computed from the sample, are removed. For further details see:
#'   \code{?boxplot.stats}
#' @param minPCT A decimal value between 0 and 1. Defines the minimum fraction
#'   of cells where the gene needs to be expressed to be included in the
#'   analysis.
#' @param maxMTratio A decimal value between 0 and 1. Defines the maximum ratio
#'   of mitochondrial reads (mitochondrial reads / library size) present in a
#'   cell to be included in the analysis. It's computed using the symbol genes
#'   starting with 'MT-' non-case sensitive.
#' @param label Optional character label prepended to progress messages when
#'   running inside a pipeline.
#' @return A dgCMatrix object with the cells and the genes that pass the quality
#'   control filters.
#' @references Ilicic, Tomislav, et al. "Classification of low quality cells from single-cell RNA-seq data." Genome biology 17.1 (2016): 29.
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
#' # Comparing dimensions before and after QC
#' dim(X)
#' dim(qcOutput)

scQC <- function(X, minLibSize = 1000, removeOutlierCells = TRUE,
                 minPCT = 0.05, maxMTratio = 0.1, label = NULL) {
  nCellsInit <- ncol(X)
  nGenesInit <- nrow(X)
  tag <- if (!is.null(label)) paste0("[", label, "] ") else ""

  # Remove negative values
  X[X < 0] <- 0

  # Filter by minimum library size
  lSize <- Matrix::colSums(X)
  X <- X[, lSize > minLibSize]
  cli::cli_alert_info(
    "{tag}Library size filter (min={minLibSize}): {ncol(X)}/{nCellsInit} cells retained"
  )

  # Remove outlier cells
  if (removeOutlierCells) {
    lSize <- Matrix::colSums(X)
    X <- X[, !lSize %in% boxplot.stats(lSize)$out]
    cli::cli_alert_info("{tag}Outlier removal: {ncol(X)} cells retained")
  }

  # Filter by mitochondrial ratio
  mtGenes <- grepl('^MT-', toupper(rownames(X)), ignore.case = TRUE)
  if (sum(mtGenes) > 0) {
    mtRate <- Matrix::colSums(X[mtGenes, ]) / Matrix::colSums(X)
    X <- X[, mtRate < maxMTratio]
    cli::cli_alert_info(
      "{tag}MT ratio filter (max={maxMTratio}): {ncol(X)} cells retained"
    )
  } else {
    cli::cli_alert_warning(
      "{tag}No mitochondrial genes found. Apoptotic cells may be present."
    )
  }

  # Filter by minimum expression percentage
  X <- X[Matrix::rowMeans(X != 0) > minPCT, ]
  cli::cli_alert_info(
    "{tag}Gene filter (minPCT={minPCT}): {nrow(X)}/{nGenesInit} genes retained"
  )

  # Convert to sparse matrix
  gNames <- rownames(X)
  cNames <- colnames(X)
  X <- as(X, 'dgCMatrix')
  rownames(X) <- gNames
  colnames(X) <- cNames

  cli::cli_alert_success("{tag}QC complete: {nrow(X)} genes x {ncol(X)} cells")
  return(X)
}
