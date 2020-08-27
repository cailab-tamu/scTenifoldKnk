library(Matrix)
manifoldAlignmentF <- function(X, Y, d = 2, lambda = 0.9){
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
  
  rm(wX,wY,wXY)
  Drow <- rowSums(W)
  Dcol <- colSums(W)
  W <- Matrix(W) 
  diag(W) <- -colSums(W) 
  
  #if(sqrt1 == 1){
  #  Dcol = Dcol^(0.5)
  #}
  
  #if(weight_method == 1){
  ## Yan's approach
  D <- Matrix::diag(Dcol)
  D <- as(D, 'dgCMatrix')
  W =  D - W * (Drow^(-1)) * Dcol
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

load('betaMALATko.RData')
set.seed(1)
MA <- manifoldAlignmentF(MALAT1$WT, MALAT1$KO, d = 2)
MA <- MA[!grepl('_Rpl|_Rps|_Rp[[:digit:]]+|_Mt-',rownames(MA), ignore.case = TRUE),]
DR <- scTenifoldKnk:::dRegulation(MA, gKO = 'Malat1')
DR$FC <- (DR$distance^2/mean(DR$distance[-c(1:2)]^2))
DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
