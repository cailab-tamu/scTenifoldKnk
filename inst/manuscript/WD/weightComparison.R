library(ggplot2)
library(igraph)
library(ggrepel)
library(statsExpressions)
library(scTenifoldNet)
library(Matrix)
plotComparison <- function(X, gene, nLabels = 15){
  diag(X$WT) <- 0
  W0 <- apply(X$WT[gene,, drop=FALSE],2,max)
  W <- W0
  set.seed(1)
  W <- W + rnorm(length(W),mean = 0, sd = 0.015)
  W <- W + rnorm(length(W),mean = 0, sd = 0.015)
  W[W0 == 0] <- rnorm(sum(W0 == 0), mean = 0, sd = 0.015)
  W[gene] <- 0
  D <- X$diffRegulation$Z
  names(D) <- X$diffRegulation$gene
  sGenes <- intersect(names(W), names(D))
  W <- (W[sGenes])
  W <- round(W,3)
  D <- D[sGenes]
  DF <- data.frame(W=W, D=D, G= sGenes)
  DF$P <- ifelse(DF$G %in% X$diffRegulation$gene[X$diffRegulation$p.value < 0.05], 8, 20)
  DF$C <- densCols(DF[,c('D','W')], colramp = hcl.colors)
  nLabels <- ifelse(nLabels < sum(X$diffRegulation$p.adj < 0.05), nLabels, sum(X$diffRegulation$p.adj < 0.05))
  DF$C[DF$G %in% X$diffRegulation$gene[X$diffRegulation$p.adj < 0.05]] <- 'red'
  DF$G[!DF$G %in% X$diffRegulation$gene[seq_len(nLabels)]] <- NA
  P <- ggplot(DF, aes(D,W, label = G)) +
    geom_point(col = DF$C, pch = DF$P) +
    geom_text_repel(segment.size = 0.2, segment.alpha = 0.5, max.iter = 1e4, force = 10, nudge_x = 0.7) +
    theme_bw() +
    xlab(Z-score~(Distance)) +
    ylab(Edge~weight) +
    labs(title = paste0(gene, collapse = ' - '), subtitle = statsExpressions::corr_test(DF, W, D, type = 'nonparametric')$expression[[1]]) +
    theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 8))
  return(P)
}

# DMD
load('../DMD/Results/GSM4116571.RData')
png('DMD.png', width = 1500, height = 1500, res = 300)
plotComparison(GSM4116571, gene = 'Dmd') + xlim(c(0,2.5))
dev.off()

# AHR
load('../AHR/Results/Preenterocytes.RData')
png('AHR.png', width = 1500, height = 1500, res = 300)
plotComparison(O, gene = 'Ahr')
dev.off()

# CFTR
load('../CFTR/Results/SRS4245406.RData')
png('CFTR.png', width = 1500, height = 1500, res = 300)
plotComparison(SRS4245406, gene = 'Cftr')
dev.off()

# MALAT1
load('../MALAT1/betaMALATko.RData')
png('MALAT1.png', width = 1500, height = 1500, res = 300)
MALAT1$diffRegulation <- scTenifoldKnk:::dRegulation(MALAT1$MA, gKO = c('A','B'))
MALAT1$diffRegulation$Z <- log1p(MALAT1$diffRegulation$Z)
plotComparison(MALAT1, gene = 'Malat1') + xlab(log~(Z-score~(Distance)+1))
dev.off()

# MECP2
png('MECP2.png', width = 1500, height = 1500, res = 300)
load('../MECP2/Results/SRS3059998.RData')
plotComparison(SRS3059998, gene = 'Mecp2')
dev.off()

# NKX21
png('NKX21.png', width = 1500, height = 1500, res = 300)
load('../NKX2-1/Results/GSM3716703.RData')
plotComparison(GSM3716703, gene = 'Nkx2-1') + xlim(c(-1.5,3.5))
dev.off()

# TREM2
png('TREM2.png', width = 1500, height = 1500, res = 300)
load('../TREM2/Results/GSE130626.RData')
plotComparison(GSE130626, gene = 'Trem2')
dev.off()

# HNF4AG
png('HNF4AG.png', width = 1500, height = 1500, res = 300)
load('../HNF4A-HNF4G/Results/GSM3477499.RData')
plotComparison(GSM3477499, gene = c('Hnf4a', 'Hnf4g')) + xlim(c(-1.5,4))
dev.off()


plotComparison <- function(X, gene, nLabels = 15){
  diag(X$WT) <- 0
  W0 <- apply(X$WT[gene,, drop=FALSE],2,max)
  W <- W0
  set.seed(1)
  W <- W + rnorm(length(W),mean = 0, sd = 0.015)
  W <- W + rnorm(length(W),mean = 0, sd = 0.015)
  W[W0 == 0] <- rnorm(sum(W0 == 0), mean = 0, sd = 0.015)
  W[gene] <- 0
  D <- X$diffRegulation$Z
  names(D) <- X$diffRegulation$gene
  sGenes <- intersect(names(W), names(D))
  W <- abs(W[sGenes])
  W <- round(W,3)
  D <- D[sGenes]
  DF <- data.frame(W=W, D=D, G= sGenes)
  DF$P <- ifelse(DF$G %in% X$diffRegulation$gene[X$diffRegulation$p.value < 0.05], 8, 20)
  DF$C <- densCols(DF[,c('D','W')], colramp = hcl.colors)
  nLabels <- ifelse(nLabels < sum(X$diffRegulation$p.adj < 0.05), nLabels, sum(X$diffRegulation$p.adj < 0.05))
  DF$C[DF$G %in% X$diffRegulation$gene[X$diffRegulation$p.adj < 0.05]] <- 'red'
  DF$G[!DF$G %in% X$diffRegulation$gene[seq_len(nLabels)]] <- NA
  P <- ggplot(DF, aes(D,W, label = G)) +
    geom_point(col = DF$C, pch = DF$P) +
    geom_text_repel(segment.size = 0.2, segment.alpha = 0.5, max.iter = 1e3, force = 15, nudge_x = 0.7) +
    theme_bw() +
    xlab(Z-score~(Distance)) +
    ylab(Abs~(Edge~weight)) +
    labs(title = paste0(gene, collapse = ' - ')) +
    theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 8))
  return(P)
}

# DMD
load('../DMD/Results/GSM4116571.RData')
png('absDMD.png', width = 1500, height = 1500, res = 300)
DMDPlot <- plotComparison(GSM4116571, gene = 'Dmd') + xlim(c(0,2.5))
print(DMDPlot)
dev.off()

# AHR
load('../AHR/Results/Preenterocytes.RData')
png('absAHR.png', width = 1500, height = 1500, res = 300)
plotComparison(O, gene = 'Ahr')
dev.off()

# CFTR
load('../CFTR/Results/SRS4245406.RData')
png('absCFTR.png', width = 1500, height = 1500, res = 300)
plotComparison(SRS4245406, gene = 'Cftr')
dev.off()

# MALAT1
load('../MALAT1/betaMALATko.RData')
png('absMALAT1.png', width = 1500, height = 1500, res = 300)
MALAT1$diffRegulation <- scTenifoldKnk:::dRegulation(MALAT1$MA, gKO = c('A','B'))
MALAT1$diffRegulation$Z <- log1p(MALAT1$diffRegulation$Z)
plotComparison(MALAT1, gene = 'Malat1') + xlab(log~(Z-score~(Distance)+1))
dev.off()

# MECP2
png('absMECP2.png', width = 1500, height = 1500, res = 300)
load('../MECP2/Results/SRS3059998.RData')
plotComparison(SRS3059998, gene = 'Mecp2')
dev.off()

# NKX21
png('absNKX21.png', width = 1500, height = 1500, res = 300)
load('../NKX2-1/Results/GSM3716703.RData')
plotComparison(GSM3716703, gene = 'Nkx2-1') + xlim(c(-1.5,3.5))
dev.off()

# TREM2
png('absTREM2.png', width = 1500, height = 1500, res = 300)
load('../TREM2/Results/GSE130626.RData')
TREM2Plot <- plotComparison(GSE130626, gene = 'Trem2')
print(TREM2Plot)
dev.off()

# HNF4AG
png('absHNF4AG.png', width = 1500, height = 1500, res = 300)
load('../HNF4A-HNF4G/Results/GSM3477499.RData')
plotComparison(GSM3477499, gene = c('Hnf4a', 'Hnf4g')) + xlim(c(-1.5,4))
dev.off()


library(patchwork)
pdf(file = 'Abs_Trem2+Dmd.pdf', width = 10, height = 5)
TREM2Plot + DMDPlot
dev.off()
