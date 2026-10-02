#' @export heatKernel
#' @title Spectral heat kernel of a gene-gene matrix
#' @description Computes the heat kernel
#'   \deqn{H = \sum_k \exp(t (\lambda_k / \lambda_{max} - 1)) v_k v_k^T}
#'   from the eigendecomposition \eqn{X = \sum_k \lambda_k v_k v_k^T} of a
#'   symmetric gene-gene matrix (a gene regulatory network or a gene-gene
#'   correlation matrix). Row \code{i} of \code{H} describes how a change in
#'   gene \code{i} diffuses over the whole network, summing paths of every
#'   length; larger \code{t} concentrates the kernel on the dominant network
#'   modules.
#' @param X A square numeric matrix (or \code{Matrix}) with genes as row and
#'   column names.
#' @param t A non-negative number. Diffusion time. \code{t = 0} returns the
#'   identity.
#' @param symmetric A boolean value (TRUE/FALSE). If TRUE, \code{X} is replaced
#'   by its symmetric part \code{(X + t(X)) / 2} before the eigendecomposition;
#'   required for directed networks such as the scTenifoldKnk WT network.
#'   Default: TRUE.
#' @return A square numeric matrix with the same dimension names as \code{X}.
#' @examples
#' set.seed(1)
#' A <- matrix(rnorm(25), 5, 5, dimnames = list(letters[1:5], letters[1:5]))
#' H <- heatKernel(A, t = 10)
#' isSymmetric(H)
heatKernel <- function(X, t = 10, symmetric = TRUE) {
  if (!is.numeric(t) || length(t) != 1 || t < 0) {
    stop("'t' must be a single non-negative number")
  }
  X <- as.matrix(X)
  if (nrow(X) != ncol(X)) {
    stop("'X' must be a square matrix")
  }
  if (isTRUE(symmetric)) {
    X <- (X + t(X)) / 2
  } else if (!isSymmetric(unname(X))) {
    stop("'X' is not symmetric; use 'symmetric = TRUE'")
  }
  E <- eigen(X, symmetric = TRUE)
  lambdaMax <- max(E$values)
  if (lambdaMax <= 0) {
    stop("The largest eigenvalue of 'X' must be positive")
  }
  w <- exp(t * (E$values / lambdaMax - 1))
  H <- E$vectors %*% (w * t(E$vectors))
  dimnames(H) <- dimnames(X)
  return(H)
}

# log1p(CPM) of the genes in 'genes' and their per-gene mean and standard deviation
.logCPM <- function(X, genes) {
  missingGenes <- setdiff(genes, rownames(X))
  if (length(missingGenes) > 0) {
    stop("The following genes are not present in 'X': ", paste(missingGenes, collapse = ", "))
  }
  libSize <- Matrix::colSums(X)
  L <- as.matrix(X[genes, , drop = FALSE])
  L <- log1p(t(t(L) / libSize) * 1e6)
  mu <- rowMeans(L)
  sdev <- sqrt(rowMeans((L - mu)^2))
  sdev[sdev == 0] <- NA
  list(L = L, mu = mu, sd = sdev)
}

# Signed change of every gene when the genes in each element of 'gKOs' are removed:
# sum over the knocked-out genes x of -mean(x) / sd(x) * H[x, g] * sd(g)
.kernelEffect <- function(H, stats, gKOs) {
  genes <- rownames(H)
  out <- t(vapply(gKOs, function(ko) {
    w <- -stats$mu[ko] / stats$sd[ko]
    w[is.na(w)] <- 0
    as.numeric(crossprod(w, H[ko, , drop = FALSE])) * stats$sd
  }, numeric(length(genes))))
  out[is.na(out)] <- 0
  dimnames(out) <- list(names(gKOs), genes)
  out
}

#' @export hkManifoldAlignment
#' @title Heat-kernel (hk) manifold alignment for in-silico knockouts
#' @description Kernel counterpart of the non-linear manifold alignment used by
#'   \code{scTenifoldKnk}. The heat kernel of the WT network is computed once
#'   (\code{\link{heatKernel}}), and the knockout of each gene (or gene set) is
#'   read from its rows: removing gene \code{x} lowers it by its mean
#'   expression and the change diffuses over the network,
#'   \deqn{\Delta_g = -\frac{mean(x)}{sd(x)} H[x, g] \, sd(g)}
#'   on log1p(CPM) expression (summed over the genes of a multi-gene knockout).
#'   The absolute value \eqn{|\Delta_g|} plays the role of the manifold
#'   alignment distance: it ranks the perturbed genes like the manifold
#'   alignment does, without recomputing an alignment for every knockout, which
#'   makes transcriptome-wide perturbation efficient.
#' @param WT The WT gene regulatory network (genes x genes, regulators as rows),
#'   as returned in \code{scTenifoldKnk(...)$tensorNetworks$WT}.
#' @param X Raw counts matrix (genes x cells) of the WT cells, with row names
#'   matching the genes of \code{WT} (typically the quality-controlled matrix
#'   used to build the network).
#' @param gKO A character vector or a list. Each element is a gene (or a
#'   character vector of genes knocked out together) to perturb. If
#'   \code{NULL}, every gene of the network is knocked out separately.
#' @param t A non-negative number. Diffusion time of the heat kernel. Default:
#'   10.
#' @param H Optional pre-computed heat kernel of \code{WT} (from
#'   \code{heatKernel(WT, t)}), to reuse it across calls.
#' @return A numeric matrix of perturbation distances \eqn{|\Delta_g|} (for
#'   the knocked-out genes themselves, their mean log1p(CPM) expression, the
#'   amount they lose) with one row per knockout (named by the knocked-out genes, joined with \code{"+"}
#'   for multi-gene knockouts) and one column per gene of the network.
#' @seealso \code{\link{heatKernel}}, \code{\link{knockoutDirection}}
hkManifoldAlignment <- function(WT, X, gKO = NULL, t = 10, H = NULL) {
  WT <- as.matrix(WT)
  genes <- rownames(WT)
  if (is.null(gKO)) gKO <- as.list(genes)
  if (!is.list(gKO)) gKO <- as.list(gKO)
  names(gKO) <- vapply(gKO, paste, character(1), collapse = "+")
  missingGenes <- setdiff(unlist(gKO), genes)
  if (length(missingGenes) > 0) {
    stop("The following genes are not present in the WT network: ", paste(missingGenes, collapse = ", "))
  }
  if (is.null(H)) {
    diag(WT) <- 0
    H <- heatKernel(WT, t = t)
  }
  S <- .logCPM(X, genes)
  D <- abs(.kernelEffect(H, S, gKO))
  # A knocked-out gene loses its whole expression: its change is its mean log1p(CPM)
  for (i in seq_along(gKO)) D[i, gKO[[i]]] <- S$mu[gKO[[i]]]
  D
}

#' @export knockoutDirection
#' @importFrom stats cor
#' @title Direction of the response to an in-silico knockout
#' @description Predicts whether each gene goes up or down after knocking out
#'   \code{gKO}, using only the WT expression data: the heat kernel
#'   (\code{\link{heatKernel}}) of the gene-gene Pearson correlation matrix of
#'   log1p(CPM) expression is diffused from the knocked-out gene(s),
#'   \deqn{s_g = -\sum_{x \in gKO} \frac{mean(x)}{sd(x)} H[x, g] \, sd(g)},
#'   and the sign of \eqn{s_g} is the predicted direction.
#'
#'   Evaluated against bulk knockdown/knockout profiles of the same cell
#'   lines, this direction is correct more often than chance, but it mostly
#'   reflects the response shared by most knockdowns along the dominant WT
#'   expression program rather than regulation specific to \code{gKO}, and its
#'   accuracy varies between cell types and WT data sets.
#' @param X Raw counts matrix (genes x cells) of the WT cells.
#' @param gKO A character vector or a list, as in
#'   \code{\link{hkManifoldAlignment}}.
#' @param genes A character vector with the genes to score (for example the
#'   genes of the WT network). Default: all genes of \code{X}.
#' @param t A non-negative number. Diffusion time of the heat kernel. Default:
#'   5.
#' @return A numeric matrix of direction scores \eqn{s_g} (positive: predicted
#'   up, negative: predicted down) with one row per knockout and one column per
#'   gene.
#' @seealso \code{\link{heatKernel}}, \code{\link{dRegulation}}
knockoutDirection <- function(X, gKO, genes = rownames(X), t = 5) {
  if (!is.list(gKO)) gKO <- as.list(gKO)
  names(gKO) <- vapply(gKO, paste, character(1), collapse = "+")
  missingGenes <- setdiff(unlist(gKO), genes)
  if (length(missingGenes) > 0) {
    stop("The following genes are not present in 'genes': ", paste(missingGenes, collapse = ", "))
  }
  S <- .logCPM(X, genes)
  R <- cor(t(S$L))
  R[is.na(R)] <- 0
  H <- heatKernel(R, t = t)
  .kernelEffect(H, S, gKO)
}
