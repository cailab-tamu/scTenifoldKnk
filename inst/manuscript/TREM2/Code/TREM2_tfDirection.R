# setwd('/data/dcosorioh/manuscript/')
#
library(Matrix)
# scTenifoldKnk <- function(countMatrix, gKO = NULL){
#   set.seed(1)
#   WT <- scTenifoldNet::makeNetworks(countMatrix, q = 0.9)
#   set.seed(1)
#   WT <- scTenifoldNet::tensorDecomposition(WT)
#   WT <- as.matrix(WT$X)
#   #KO <- rCUR::CUR(WT, sv = RSpectra::svds(WT, 5))
#   #C <- KO@C
#   #C[gKO,] <- 0
#   #KO <- C %*% KO@U %*% KO@R
#   KO <- WT
#   KO[gKO,] <- 0
#   set.seed(1)
#   MA <- scTenifoldNet::manifoldAlignment(WT, KO)
#   set.seed(1)
#   DR <- scTenifoldNet::dRegulation(MA)
#   outputList <- list()
#   outputList$WT <- WT
#   outputList$KO <- KO
#   outputList$diffRegulation <- DR
#   return(outputList)
# }
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')

WT <- readMM('../Data/matrix_WT.mtx')
rownames(WT) <- readLines('../Data/geneNames_WT.txt')
colnames(WT) <- readLines('../Data/barcodes_WT.txt')
WT <- scQC(WT)
WT <- WT[rowMeans(WT != 0) > 0.05,]
WT <- WT[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(WT), ignore.case = TRUE),]
rownames(WT) <- toupper(rownames(WT))
# WT <- scTenifoldKnk(WT, gKO = 'Trem2')
# save(WT, file = 'TREM2_GSE130626.RData')

library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
library(scTenifoldNet)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa2mmu_SYMBOL.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/graphs/strictDirection.R')


WT <- WT[unique(unlist(ENCODE)),]
#W <- pcNet(as.matrix(WT))
load('W.RData')


ENCODE <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=ENCODE_TF_ChIP-seq_2015')
geneList <- toupper(rownames(W))
tfList <- lapply(strsplit(names(ENCODE), split = ' '), function(X){X[1]})
tfList <- unlist(tfList)
ENCODE <- ENCODE[tfList %in% geneList]
tfList <- tfList[tfList %in% geneList]
ENCODE <- lapply(ENCODE, function(X){X[X %in% geneList]})
tfList <- tfList[order(lengths(ENCODE))]
ENCODE <- ENCODE[order(lengths(ENCODE))]

# ENCODE <- lapply(unique(tfList), function(TF){
#   unique(unlist(ENCODE[tfList %in% TF]))
# })
# tfList <- unique(tfList)


sapply(c(0, 0.25, 0.5, 0.75, 1), function(lambda){
  W <- strictDirection(W, lambda = lambda)
  library(pbapply)
  dW <- pblapply(seq_along(tfList), function(X){
    I <- W[ENCODE[[X]],tfList[X]]
    O <- W[tfList[X],ENCODE[[X]]]
    #mean(O - I)
    I <- data.frame(TF = tfList[X], D = 'I', W = IQR(I))
    O <- data.frame(TF = tfList[X], D = 'O', W = IQR(O))
    R <- rbind(I,O)
    return(R)
  })
  dW <- do.call(rbind.data.frame, dW)
  dW <- dW[order(dW$W, decreasing = TRUE),]
  write.csv(dW, file = paste0('tf-targetALLlambda',lambda,'.csv'))
  T1 <- ks.test(dW$W[dW$D == 'O'], dW$W[dW$D == 'I'], alternative = 'g', exact = TRUE)
  T2 <- ks.test(dW$W[dW$D == 'I'], dW$W[dW$D == 'O'], alternative = 'g', exact = TRUE)
  library(ggplot2)
  png(paste0('tf-targetALLlambda',lambda,'.png'), width = 2000, height = 1500, res = 300)
  print(ggplot(dW, aes(W, color = D)) +
    stat_ecdf() + theme_bw() +
    ylab('Empirical Cumulative Density Function (ECDF)') +
    xlab('IQR(TF-Targets Weight)') +
    labs(color = 'Direction', title = parse(text = paste0('lambda==',lambda)), subtitle = parse(text = paste0('One-side~K-S~test~O>I~P-value == ',formatC(T2$p.value, digits = 3)))))
  dev.off()
  png(paste0('tf-targetTOP500lambda',lambda,'.png'), width = 2000, height = 1500, res = 300)
  print(ggplot(dW[1:500,], aes(W, color = D)) +
          stat_ecdf() + theme_bw() +
          ylab('Empirical Cumulative Density Function (ECDF)') +
          xlab('IQR(TF-Targets Weight)') +
          labs(color = 'Direction', title = parse(text = paste0('lambda==',lambda)), subtitle = parse(text = paste0('One-side~K-S~test~O>I~P-value == ',formatC(T2$p.value, digits = 3)))))
  dev.off()
  png(paste0('tf-targetBPlambda',lambda,'.png'), width = 1000, height = 1500, res = 300)
  print(ggplot(dW, aes(D, W, color = D)) + geom_boxplot(outlier.colour = NA) + theme_bw() + xlab('Direction') + ylab('IQR(TF-Targets Weight)') + geom_jitter(alpha = 0.1) + labs(color = 'Direction'))
  dev.off()
})

