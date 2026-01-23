qFilter <- function(X, q = 0) {
  X <- as.matrix(X)
  X[abs(X) < stats::quantile(abs(X), 0.95)] <- 0
  X <- Matrix::Matrix(X)
  return(X)
}
