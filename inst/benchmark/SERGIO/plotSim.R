setwd("~/scTenifoldKNK/inst/benchmark/SERGIO")
A <- read.csv('A.csv', header = FALSE)
rownames(A) <- paste0('g',seq_len(nrow(A)))
colnames(A) <- paste0('C', seq_len(ncol(A)))
# Q <- 0.8
# gKO <- c(1)
scTenifoldKnk <- function(countMatrix, gKO = NULL){
  set.seed(1)
  WT <- scTenifoldNet::makeNetworks(countMatrix, q = Q)
  set.seed(1)
  WT <- scTenifoldNet::tensorDecomposition(WT)
  WT <- as.matrix(WT$X)
  #KO <- rCUR::CUR(WT, sv = RSpectra::svds(WT, 5))
  #C <- KO@C
  #C[gKO,] <- 0
  #KO <- C %*% KO@U %*% KO@R
  KO <- WT
  KO[gKO,] <- 0
  set.seed(1)
  MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2)
  set.seed(1)
  DR <- scTenifoldNet::dRegulation(MA)
  D <- DR$distance[-seq_along(gKO)]
  DR$FC <- DR$distance^2/mean(D^2)
  DR$p.value <- c(rep(0,length(gKO)),pchisq((D^2/mean(D^2)), df = 1, lower.tail = FALSE))
  DR$p.adj <- p.adjust(DR$p.value)
  outputList <- list()
  outputList$WT <- Matrix::Matrix(WT)
  outputList$KO <- Matrix::Matrix(KO)
  outputList$manifoldAlignment <- MA
  outputList$diffRegulation <- DR
  return(outputList)
}

Q <- 0.9
K1 <- scTenifoldKnk(A,1)
K1$diffRegulation <- K1$diffRegulation[paste0(1:10),]
K7 <- scTenifoldKnk(A,7)
K7$diffRegulation <- K7$diffRegulation[paste0(1:10),]

plot(K1$diffRegulation$FC)
g <- c(rep(1,6),rep(2,4))
par(mfrow=c(1,2))
boxplot(K1$diffRegulation$Z ~ g)
boxplot(K7$diffRegulation$Z ~ g)
