#' @import Matrix
manifoldAlignment <- function(X, Y, d = 2, lambda = 0.9){
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  n = dim(X)[1]
  Y <- Y[sharedGenes, sharedGenes]
  L <- diag(length(sharedGenes))
  wX <- X+1
  wY <- Y+1
  wXY <- lambda * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  diag(W) <- 0
  
  Drow <- apply(W, 1, sum)
  Dcol <- apply(W, 2, sum)
  W <- Matrix(W)  
  
  #if(sqrt1 == 1){
  #  Dcol = Dcol^(0.5)
  #}
  
  #if(weight_method == 1){
  ## Yan's approach
  W = diag(Dcol) - W * (Drow^(-1)) * Dcol
  W = t(W) %*% W
  #} else {
  #  Dcol[c(1:n)+n] = Dcol[c(1:n)]
  #  W = diag(Dcol) - W * (Drow^(-1)) * Dcol
  #}
  
 
  E <- suppressWarnings(RSpectra::eigs(W, d*2, 'SR'))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[,newOrder]
  E$vectors <- E$vectors[,E$values > 1e-8]
  alignedNet <- E$vectors[,seq_len(d)]
  colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}


#' @import Matrix test by Yan
#' @ k: 2,3,4, gamma: 0~1. gamma = 1 means using the same weight for all order. methoduse = 1, use cumulant orders, not 1, only kth order.
manifoldAlignment_test <- function(X, Y, d = 2, lambda = 0.9, k = 2, gamma = 0.5, methoduse = 1){
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  n = dim(X)[1]
  Y <- Y[sharedGenes, sharedGenes]
  
  Xuse = X
  Yuse = Y
  Xho = X
  Yho = Y
  if(k >=2){
    for(i in 1:(k-1)){
      Xho = Xho %*% X
      Yho = Yho %*% Y
      Xuse = Xuse + Xho * gamma^i
      Yuse = Yuse + Yho * gamma^i
    }
  }
  
  if(methoduse == 1){
    X = Xuse
    Y = Yuse
  }else{
    X = Xho
    Y = Yho
  }
  
  L <- diag(length(sharedGenes))
  wX <- X+1
  wY <- Y+1
  wXY <- lambda * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  diag(W) <- 0
  
  Drow <- apply(W, 1, sum)
  Dcol <- apply(W, 2, sum)
  W <- Matrix(W)  
  
  #if(sqrt1 == 1){
  #  Dcol = Dcol^(0.5)
  #}
  
  #if(weight_method == 1){
  ## Yan's approach
  W = diag(Dcol) - W * (Drow^(-1)) * Dcol
  W = t(W) %*% W
  #} else {
  #  Dcol[c(1:n)+n] = Dcol[c(1:n)]
  #  W = diag(Dcol) - W * (Drow^(-1)) * Dcol
  #}
  
 
  E <- suppressWarnings(RSpectra::eigs(W, d*2, 'SR'))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[,newOrder]
  E$vectors <- E$vectors[,E$values > 1e-8]
  alignedNet <- E$vectors[,seq_len(d)]
  colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}
