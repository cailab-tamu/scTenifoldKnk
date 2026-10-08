#' @title Map of transcriptome-wide perturbation profiles
#' @description Places every virtual knockout from a transcriptome-wide run of
#'   \code{\link{scTenifoldKnk}} on a two-dimensional map together with a
#'   reference signature, such as a disease-vs-control expression signature, and
#'   highlights the genes whose knockout reproduces it. Each knockout is
#'   described by its signed perturbation profile: the rank of the perturbation
#'   distance of every gene, scaled to (0, 1], multiplied by its predicted
#'   direction. The similarity of a profile to the signature is the cosine
#'   between both over the genes they share; the knocked-out gene itself is left
#'   out of its own profile.
#' @param X A list. Output from \code{\link{scTenifoldKnk}} with
#'   \code{transcriptomeWide = TRUE} and \code{dr_direction = TRUE}, which
#'   contains \code{perturbationDistances} and \code{perturbationDirections}.
#' @param signature A named numeric vector, for example the log fold changes of
#'   a disease-vs-control comparison in the same cell type. Names must be gene
#'   identifiers matching those of the network.
#' @param genes Character. Optional genes to highlight, such as known causal
#'   genes. Default: \code{NULL}.
#' @param layout Character, \code{"signature"} or \code{"pca"}.
#'   \code{"signature"} (default) places each knockout by its similarity to the
#'   signature (x axis) and by the main axis of the remaining variation of the
#'   profiles (y axis). \code{"pca"} shows the first two principal components of
#'   the profiles and projects the signature into the same space.
#' @param nLabels Integer. Number of knockouts with the highest similarity to
#'   label, in addition to \code{genes}. Default: 10.
#' @param plot Logical. If \code{TRUE} (default), draw the map.
#' @return Invisibly, a data frame with one row per knockout, sorted by
#'   decreasing similarity: \code{gene}, \code{similarity} (cosine to the
#'   signature), \code{rank} (1 = most similar), \code{percentile} (fraction of
#'   knockouts at least as similar), and the map coordinates \code{x} and
#'   \code{y}. With \code{layout = "pca"}, the coordinates of the signature are
#'   stored in the \code{"signature"} attribute.
#' @examples
#' \donttest{
#' library(scTenifoldKnk)
#'
#' # Counts with three gene modules
#' set.seed(1)
#' modules <- matrix(rgamma(3 * 600, shape = 2), 3, 600)
#' loadings <- matrix(0, 60, 3)
#' loadings[cbind(1:60, rep(1:3, length.out = 60))] <- runif(60, 0.5, 2)
#' X <- matrix(rpois(60 * 600, lambda = 5 * loadings %*% modules), 60, 600)
#' rownames(X) <- paste0('g', 1:60)
#' colnames(X) <- paste0('c', 1:600)
#'
#' # Knock out every gene of the network
#' O <- scTenifoldKnk(X, transcriptomeWide = TRUE, qc = FALSE,
#'                    nc_nNet = 3, nc_nCells = 300, nCores = 1)
#'
#' # A signature that resembles the knockout of g1
#' sig <- O$perturbationDistances['g1', ] * O$perturbationDirections['g1', ]
#' sig <- sig + rnorm(length(sig), sd = sd(sig))
#'
#' M <- perturbationMap(O, signature = sig, genes = 'g1')
#' head(M)
#' }
#' @export
#' @importFrom graphics abline points text
#' @importFrom grDevices hcl.colors
#' @importFrom graphics legend
#' @importFrom stats quantile
perturbationMap <- function(X, signature, genes = NULL, layout = c("signature", "pca"),
                            nLabels = 10, plot = TRUE) {
  layout <- match.arg(layout)
  D <- X$perturbationDistances
  S <- X$perturbationDirections
  if (is.null(D) || is.null(S)) {
    stop("X must be the output of scTenifoldKnk(transcriptomeWide = TRUE, dr_direction = TRUE)")
  }
  if (!is.numeric(signature) || is.null(names(signature))) {
    stop("signature must be a named numeric vector")
  }
  D <- as.matrix(D)
  D[is.na(D)] <- 0

  # Signed perturbation profiles: rank of the distance (0, 1] times the predicted direction
  P <- t(apply(D, 1, rank)) / ncol(D) * as.matrix(S)[, colnames(D), drop = FALSE]
  dimnames(P) <- dimnames(D)
  self <- match(rownames(P), colnames(P))
  P[cbind(which(!is.na(self)), self[!is.na(self)])] <- 0

  shared <- intersect(colnames(P), names(signature)[is.finite(signature)])
  if (length(shared) < 10) {
    stop("signature shares fewer than 10 genes with the network")
  }
  P <- P[, shared, drop = FALSE]
  s <- signature[shared]
  nP <- sqrt(rowSums(P^2))
  P <- P[nP > 0, , drop = FALSE]
  nP <- nP[nP > 0]
  if (nrow(P) < 3) {
    stop("at least 3 knockouts with a non-zero profile are needed")
  }
  similarity <- drop(P %*% s) / (nP * sqrt(sum(s^2)))

  if (layout == "signature") {
    # x: similarity to the signature; y: main axis of the profiles orthogonal to the signature
    u <- s / sqrt(sum(s^2))
    R <- P / nP
    R <- R - outer(drop(R %*% u), u)
    R <- sweep(R, 2, colMeans(R))
    y <- svd(R, nu = 1, nv = 0)$u[, 1]
    xy <- cbind(similarity, y)
  } else {
    # First two principal components of the profiles; the signature is projected and scaled
    # to the extent of the knockouts
    center <- colMeans(P)
    sv <- svd(sweep(P, 2, center), nu = 2, nv = 2)
    xy <- sv$u %*% diag(sv$d[1:2], 2)
    sigXY <- drop((s - center) %*% sv$v)
    sigXY <- sigXY / sqrt(sum(sigXY^2)) * quantile(sqrt(rowSums(xy^2)), 0.99)
  }

  out <- data.frame(gene = rownames(P), similarity = similarity,
                    rank = rank(-similarity, ties.method = "min"),
                    percentile = rank(-similarity, ties.method = "max") / length(similarity),
                    x = xy[, 1], y = xy[, 2], row.names = NULL)
  if (layout == "pca") attr(out, "signature") <- c(x = sigXY[[1]], y = sigXY[[2]])

  missingGenes <- setdiff(genes, out$gene)
  if (length(missingGenes) > 0) {
    warning("The following genes are not among the perturbed genes with a non-zero profile: ",
            paste(missingGenes, collapse = ", "))
  }

  if (isTRUE(plot)) {
    lim <- max(abs(out$similarity))
    pal <- hcl.colors(101, "Blue-Red 3")
    col <- pal[round((out$similarity / lim + 1) * 50) + 1]
    xr <- range(c(out$x, if (layout == "pca") sigXY[1]))
    yr <- range(c(out$y, if (layout == "pca") sigXY[2]))
    plot(out$x, out$y, pch = 16, cex = 0.6, col = col, xlim = xr, ylim = yr,
         xlab = if (layout == "signature") "Similarity to the signature (cosine)" else "PC1",
         ylab = if (layout == "signature") "Main axis of the remaining variation" else "PC2",
         main = "Perturbation map")
    if (layout == "signature") abline(v = 0, lty = 2, col = "grey60")
    top <- setdiff(out$gene[order(-out$similarity)][seq_len(min(nLabels, nrow(out)))], genes)
    i <- match(top, out$gene)
    j <- match(intersect(genes, out$gene), out$gene)
    if (length(j) > 0) points(out$x[j], out$y[j], pch = 1, cex = 1.6, lwd = 1.5)
    if (layout == "pca") points(sigXY[1], sigXY[2], pch = 23, cex = 2, bg = "gold")

    # legend in the corner with the fewest knockouts, avoiding the labelled ones
    left <- out$x < mean(xr)
    low <- out$y < mean(yr)
    w <- 1 + 100 * seq_len(nrow(out)) %in% c(i, j)
    corners <- c(topleft = sum(w[left & !low]), topright = sum(w[!left & !low]),
                 bottomleft = sum(w[left & low]), bottomright = sum(w[!left & low]))
    leg <- legend(names(which.min(corners)), bty = "n", cex = 0.7, pch = c(16, 16, 1),
                  pt.cex = c(0.8, 0.8, 1.3), col = c(pal[101], pal[1], "black"),
                  legend = c("Reproduces the signature", "Opposes the signature",
                             if (length(j) > 0) "Highlighted genes" else NA)[c(TRUE, TRUE, length(j) > 0)])$rect

    # labels: highlighted genes first, then the most similar knockouts, then the signature
    labX <- c(out$x[j], out$x[i], if (layout == "pca") sigXY[1])
    labY <- c(out$y[j], out$y[i], if (layout == "pca") sigXY[2])
    labs <- c(out$gene[j], out$gene[i], if (layout == "pca") "signature")
    labCex <- c(rep(0.8, length(j)), rep(0.6, length(i)), if (layout == "pca") 0.8)
    labFont <- c(rep(2, length(j)), rep(1, length(i)), if (layout == "pca") 2)
    labCol <- c(rep("black", length(j)), rep("grey30", length(i)), if (layout == "pca") "black")
    labGap <- c(rep(2, length(j)), rep(1, length(i)), if (layout == "pca") 2.5)
    placeLabels(labX, labY, labs, cex = labCex, font = labFont, col = labCol, gap = labGap,
                avoid = rbind(c(leg$left, leg$top - leg$h, leg$left + leg$w, leg$top)))
  }
  invisible(out[order(-out$similarity), , drop = FALSE])
}

# Place text labels next to their points without overlapping each other, the labelled points or the
# boxes in 'avoid' (rows: x0, y0, x1, y1 in user coordinates); 'gap' scales the clearance around each
# point, for markers larger than a plain point. Candidate positions are tried around each point at
# increasing distances; the first one free of overlaps is used, otherwise the least overlapping.
# Labels placed away from their point are joined to it by a thin line.
#' @importFrom graphics par segments strheight strwidth text
#' @noRd
placeLabels <- function(x, y, labels, cex, font, col, gap = rep(1, length(labels)), avoid = NULL) {
  if (length(labels) == 0) return(invisible(NULL))
  usr <- par("usr")
  cw <- strwidth("M", cex = 0.7)
  ch <- strheight("M", cex = 0.7)
  placed <- rbind(avoid, cbind(x - gap * cw / 2, y - gap * ch / 2, x + gap * cw / 2, y + gap * ch / 2))
  overlap <- function(b) {
    sum(pmax(0, pmin(b[3], placed[, 3]) - pmax(b[1], placed[, 1])) *
          pmax(0, pmin(b[4], placed[, 4]) - pmax(b[2], placed[, 2])))
  }
  dirs <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1), c(1, 1), c(-1, 1), c(1, -1), c(-1, -1))
  for (k in seq_along(labels)) {
    w <- strwidth(labels[k], cex = cex[k], font = font[k])
    h <- strheight(labels[k], cex = cex[k], font = font[k]) * 1.3
    best <- NULL
    bestScore <- Inf
    for (step in c(1, 2, 3.5, 5, 7)) {
      for (d in seq_len(nrow(dirs))) {
        cx <- x[k] + dirs[d, 1] * (w / 2 + cw * 0.6 * step * gap[k])
        cy <- y[k] + dirs[d, 2] * (h / 2 + ch * 0.6 * step * gap[k])
        b <- c(cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
        outside <- b[1] < usr[1] || b[3] > usr[2] || b[2] < usr[3] || b[4] > usr[4]
        score <- overlap(b) / (w * h) + 10 * outside
        if (score < bestScore) {
          bestScore <- score
          best <- c(cx, cy, step, dirs[d, ])
        }
      }
      if (bestScore == 0) break
    }
    if (best[3] > 1) {
      segments(x[k], y[k], best[1] - best[4] * w / 2, best[2] - best[5] * h / 2,
               col = "grey60", lwd = 0.6)
    }
    text(best[1], best[2], labels[k], cex = cex[k], font = font[k], col = col[k])
    placed <- rbind(placed, c(best[1] - w / 2, best[2] - h / 2, best[1] + w / 2, best[2] + h / 2))
  }
  invisible(NULL)
}
