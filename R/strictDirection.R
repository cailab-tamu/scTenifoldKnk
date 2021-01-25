strictDirection <- function(X, lambda = 1){
  S <- as.matrix(X)
  S[abs(S) < abs(t(S))] <- 0
  O <- (((1-lambda) * X) + (lambda * S))
  O <- Matrix::Matrix(O)
  return(O)
}
