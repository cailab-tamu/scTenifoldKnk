#' @importFrom stats median
#' @importFrom Matrix Matrix
vstNorm <- function(X){
  X <- as.matrix(X)
  X[X == 0] <- NA
  loggeomeans <- rowMeans(log(X), na.rm = TRUE)
  sf <- apply(X, 2, function(cnts){
    exp(median((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0], na.rm = TRUE))
  })
  X <- t(t(X)/sf)
  X[is.na(X)] <- 0
  X <- Matrix(X)
  return(X)
}
