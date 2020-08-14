library(ggplot2)
library(statsExpressions)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(pbapply)
library(Matrix)
plotDR <- function(X,
                   labelGenes = 'FDR',
                   boldGenes = NULL,
                   title = NULL,
                   subtitle = NULL) {
  require(ggplot2)
  require(ggrepel)
  # if(!any(X$diffRegulation$distance > 1e-10)){
  #   X$diffRegulation$p.value[X$diffRegulation$distance <= 1e-10] <- 1
  #   X$diffRegulation$p.adj[X$diffRegulation$distance <= 1e-10] <- 1  
  # }
  o <- -log10(X$diffRegulation$p.value)
  e <- -log10(pchisq(
    sort(rchisq(nrow(X$diffRegulation), df = 1), decreasing = TRUE),
    df = 1,
    lower.tail = FALSE
  ))
  dF <- data.frame(X = e,
                   Y = o,
                   geneID = X$diffRegulation$gene)
  dF <- as.data.frame.array(dF)
  if ('FDR' %in% labelGenes[1]) {
    dF$geneID[X$diffRegulation$p.adj > 0.05] = ''
  }
  if ('P' %in% labelGenes[1]) {
    dF$geneID[X$diffRegulation$p.value > 0.05] = ''
  }
  if (length(labelGenes) > 1) {
    dF$geneID[!X$diffRegulation$gene %in% labelGenes] = ''
  }
  geneColor <-
    ifelse(X$diffRegulation$p.value < 0.05, 'red', 'black')
  genePoint <- 16
  plotQQ <- ggplot(dF, aes(X, Y, label = geneID)) +
    geom_point(color = geneColor, pch = genePoint) +
    theme_bw() +
    geom_text_repel(
      segment.color = 'gray60',
      segment.alpha = 0.5,
      max.iter = 1e4,
      aes(fontface = ifelse(dF$geneID %in% boldGenes, 2, 1)),
      box.padding = .08
    ) +
    labs(
      y = expression(-log[1 * 0] * " (Observed P-values)"),
      x = expression(-log[1 * 0] * " (Expected P-values)"),
      title = title,
      subtitle = subtitle
    )
  return(plotQQ)
}

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

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
plotKO <- function(gKO){
  Y <- X
  Y[gKO,] <- 0
  MA <- scTenifoldNet::manifoldAlignment(X,Y, d = 2)
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
png('fig2_BenchmarkA.png', width = 2000, height = 2500, res = 300)
A/B/C
dev.off()


plotKO <- function(gKO){
  X <- countMatrix
  Y <- countMatrix
  Y[gKO,] <- 0
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  set.seed(1)
  X <- scTenifoldNet::makeNetworks(X, nNet = 20, q = 0.8)
  set.seed(1)
  Y <- scTenifoldNet::makeNetworks(Y, nNet = 20, q = 0.8)
  set.seed(1)
  X <- scTenifoldNet::tensorDecomposition(X)
  set.seed(1)
  Y <- scTenifoldNet::tensorDecomposition(Y)
  set.seed(1)
  MA <- scTenifoldNet::manifoldAlignment(X$X,Y$X, d = 30)
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
png('fig2_BenchmarkX.png', width = 2000, height = 2500, res = 300)
A/B/C
dev.off()