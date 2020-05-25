library(ggplot2)
library(statsExpressions)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(pbapply)
manifoldAlignment_yan <- function(X, Y, d = 30, weight_method = 1, lambda = 0.9, sqrt1 = 0){
  library(Matrix)
  library(RSpectra)
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
  W <- as(W, "sparseMatrix")  
  
  if(sqrt1 == 1){
    Dcol = Dcol^(0.5)
  }
  
  if(weight_method == 1){
    W = diag(Dcol) - W * (Drow^(-1)) * Dcol
    #  diag(Dcol) -diag(Dcol) %*% diag(Drow^(-1)) %*% W
  }else{
    Dcol[c(1:n)+n] = Dcol[c(1:n)]
    W = diag(Dcol) - W * (Drow^(-1)) * Dcol
  }
  
  W = t(W) %*% W
  
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
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')
source('~/../../Downloads/new_manifoldalignment.R')

SERGIO <- read.csv('SERGIO/simulationOutput.csv', header = FALSE)
rownames(SERGIO) <- paste0('G', seq_len(nrow(SERGIO)))
colnames(SERGIO) <- paste0('C', seq_len(ncol(SERGIO)))

mean(SERGIO != 0)
corValues <- as.matrix(scTenifoldNet::pcNet(as.matrix(SERGIO), symmetric = TRUE, nComp = 5))
png('fig1_PCR.png', width = 2300, height = 2000, res = 300)
ComplexHeatmap::Heatmap(corValues, row_order = seq_len(nrow(SERGIO)), column_order = seq_len(nrow(SERGIO)), name = 'PCR', show_row_names = FALSE, show_column_names = FALSE) +
  rowAnnotation(link = anno_mark(at = c(20,50,100),
                                 labels = paste0('G',c(20,50,100))))
dev.off()

countMatrix <- SERGIO
set.seed(1)
X <- scTenifoldNet::makeNetworks(countMatrix, q = 0.9)
X <- scTenifoldNet::tensorDecomposition(X)
X <- X$X
X <- (X + t(X))/2
plotKO <- function(gKO){
  Y <- X
  Y[gKO,] <- 0
  MA <- manifoldAlignment_yan(X,Y, d = 2)
  DR <- scTenifoldNet::dRegulation(MA)
  #DR$distance[DR$distance < 1e-12] <- 1e-12
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[paste0(1:100),]
  DR$C <- paste0('C', c(rep(1,5),rep(2,10),rep(3,25), rep(4,40), rep(5,20)))
  TPR <- round(mean(DR$p.value[DR$C %in% DR[gKO,]$C] < 0.05),2)
  TNR <- round(1-mean(DR$p.value[!DR$C %in% DR[gKO,]$C] < 0.05),2)
  BA <- round(mean(c(TPR,TNR)),2)
  A <- ggplot(DR, aes(x = 1:100, y = Z)) +
    geom_point(color = ifelse(DR$p.value < 0.05, 'red', 'black')) +
    theme_minimal() + xlab('Gene') + ylab('Z-score(Distance)') +
    geom_vline(xintercept = 5.5, lty = 2, col = 'gray50') +
    geom_vline(xintercept = 15.5, lty = 2, col = 'gray50') +
    geom_vline(xintercept = 40.5, lty = 2, col = 'gray50') +
    geom_vline(xintercept = 80.5, lty = 2, col = 'gray50') +
    labs(title = paste0('G', gKO, ' Knockout'),
         subtitle = paste0('Sensitivity = ', TPR,' Specificity = ', TNR,'\nBalanced Accuracy = ', BA))
  DR <- DR[order(DR$distance, decreasing = TRUE),]
  B <- plotDR(list(diffRegulation = DR), labelGenes = 'P') + theme_minimal()
  return(A + B)
}
A <- plotKO(20)
B <- plotKO(50)
C <- plotKO(100)
png('fig2_BenchmarkYan.png', width = 2000, height = 2500, res = 300)
A/B/C
dev.off()
