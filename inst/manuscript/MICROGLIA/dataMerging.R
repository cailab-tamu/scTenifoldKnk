library(MASS)
library(ComplexHeatmap)
library(pbapply)
library(Seurat)
library(Matrix)
library(ggplot2)
library(RSpectra)
library(Rtsne)
library(igraph)
library(fgsea)
library(patchwork)
library(ggrepel)
library(ggalt)

zKO <- read.csv('dataMerged.csv', row.names = 1)
zKO <- as.matrix(zKO)
zKO <- zKO[,complete.cases(t(zKO))]
sGenes <- intersect(rownames(zKO), colnames(zKO))
zKO <- zKO[sGenes, sGenes]
zKO <- zKO[,diag(zKO) > 0]

S <- cov(zKO)
S <- eigs(S, 50)
S <- S$vectors
set.seed(1)
S <- Rtsne(S)$Y
rownames(S) <- colnames(zKO)

eMatrix <- Read10X('../filtered_gene_bc_matrices_mex/mm10')
annotations_eMatrix <- read.csv('../annot_JP34353637.csv', stringsAsFactors = FALSE)
eMatrix <- eMatrix[,annotations_eMatrix$cell[annotations_eMatrix$cluster %in% 'WT microglia']]

eMatrix <- (t(t(eMatrix)/colSums(eMatrix)))* 1e4
eMatrix <- t(eMatrix)
eMatrix <- eMatrix[,colnames(zKO)]
eMatrix <- cov(as.matrix(eMatrix))
eMatrix <- eigs(eMatrix, 50)
eMatrix <- eMatrix$vectors
set.seed(1)
eMatrix <- Rtsne(eMatrix)$Y
rownames(eMatrix) <- colnames(zKO)

aMatrix <- read.csv('../matrixA0.csv', header=FALSE)
geneList <- readLines('../genesA0.csv')[-1]
rownames(aMatrix) <- colnames(aMatrix) <- geneList
aMatrix <- aMatrix[colnames(zKO), colnames(zKO)]
aMatrix <- cov(as.matrix(aMatrix))
aMatrix <- eigen(aMatrix, 50)
aMatrix <- aMatrix$vectors
set.seed(1)
aMatrix <- Rtsne(aMatrix)$Y
rownames(aMatrix) <- colnames(zKO)


KEGG <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
GOBP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018')
BIOP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')


makeECDF <- function(db, title){
  tsnePositions <- S#@reductions$tsne@cell.embeddings
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  avgDistanceKNK <- lapply(db, function(gSet){
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    if(length(gSet) > 2){
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      mean(gSet)  
    }
  })
  
  tsnePositions <- eMatrix#@reductions$tsne@cell.embeddings
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  avgDistanceEXP <- lapply(db, function(gSet){
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    if(length(gSet) > 2){
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      mean(gSet)  
    }
  })
  tsnePositions <- aMatrix#@reductions$tsne@cell.embeddings
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  avgDistanceNET <- lapply(db, function(gSet){
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    if(length(gSet) > 2){
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      mean(gSet)  
    }
  })
  
  gSets <- intersect(names(avgDistanceEXP), names(avgDistanceKNK))
  gSets <- intersect(gSets, names(avgDistanceNET))
  
  
  expDist <- unlist(avgDistanceEXP[gSets])
  knkDist <- unlist(avgDistanceKNK[gSets])
  netDist <- unlist(avgDistanceNET[gSets])
  
  expDist <- data.frame(dist = expDist, method = 'Expression')
  knkDist <- data.frame(dist = knkDist, method = 'scTenifoldKnk')
  netDist <- data.frame(dist = netDist, method = 'Network')
  allDist <- rbind(expDist, knkDist, netDist)
  allDist$method <- factor(allDist$method, levels = c('Expression', 'Network', 'scTenifoldKnk'))
  
  kTest <- kruskal.test(allDist$dist, allDist$method)
  kTest <- formatC(kTest$p.value, digits = 3)
  
  ggplot(allDist, aes(dist, color = method)) + stat_ecdf() + 
    theme_bw() + 
    xlab(Avg~(Mahalinobis~Distance)) + 
    ylab('Empirical Cumulative\nDensity Function (ECDF)') + 
    labs(title = title, subtitle = paste0('Kruskal test P-value = ', kTest)) + 
    theme(plot.title = element_text(face = 2, size = 25)) +
    guides(color=guide_legend(title="Method:"))
}
makeiqrECDF <- function(db, title){
  tsnePositions <- S#@reductions$tsne@cell.embeddings
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  avgDistanceKNK <- lapply(db, function(gSet){
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    if(length(gSet) > 2){
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      IQR(gSet)  
    }
  })
  
  tsnePositions <- eMatrix#@reductions$tsne@cell.embeddings
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  avgDistanceEXP <- lapply(db, function(gSet){
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    if(length(gSet) > 2){
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      IQR(gSet)  
    }
  })
  tsnePositions <- aMatrix#@reductions$tsne@cell.embeddings
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  avgDistanceNET <- lapply(db, function(gSet){
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    if(length(gSet) > 2){
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      IQR(gSet)  
    }
  })
  
  gSets <- intersect(names(avgDistanceEXP), names(avgDistanceKNK))
  gSets <- intersect(gSets, names(avgDistanceNET))
  
  
  expDist <- unlist(avgDistanceEXP[gSets])
  knkDist <- unlist(avgDistanceKNK[gSets])
  netDist <- unlist(avgDistanceNET[gSets])
  
  expDist <- data.frame(dist = expDist, method = 'Expression')
  knkDist <- data.frame(dist = knkDist, method = 'scTenifoldKnk')
  netDist <- data.frame(dist = netDist, method = 'Network')
  allDist <- rbind(expDist, knkDist, netDist)
  allDist$method <- factor(allDist$method, levels = c('Expression', 'Network', 'scTenifoldKnk'))
  
  kTest <- kruskal.test(allDist$dist, allDist$method)
  kTest <- formatC(kTest$p.value, digits = 3)
  
  ggplot(allDist, aes(dist, color = method)) + stat_ecdf() + 
    theme_bw() + 
    xlab(IQR~(Mahalinobis~Distance)) + 
    ylab('Empirical Cumulative\nDensity Function (ECDF)') + 
    labs(title = title, subtitle = paste0('Kruskall test P-value = ', kTest)) + 
    theme(plot.title = element_text(face = 2, size = 25)) +
    guides(color=guide_legend(title="Method:"))
}
makeiqrBP <- function(db, title){
  tsnePositions <- S#@reductions$tsne@cell.embeddings
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  avgDistanceKNK <- lapply(db, function(gSet){
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    if(length(gSet) > 2){
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      IQR(gSet)  
    }
  })
  
  tsnePositions <- eMatrix#@reductions$tsne@cell.embeddings
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  avgDistanceEXP <- lapply(db, function(gSet){
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    if(length(gSet) > 2){
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      IQR(gSet)  
    }
  })
  
  tsnePositions <- aMatrix#@reductions$tsne@cell.embeddings
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  avgDistanceNET <- lapply(db, function(gSet){
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    if(length(gSet) > 2){
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      IQR(gSet)  
    }
  })
  
  gSets <- intersect(names(avgDistanceEXP), names(avgDistanceKNK))
  gSets <- intersect(gSets, names(avgDistanceNET))
  
  
  expDist <- unlist(avgDistanceEXP[gSets])
  knkDist <- unlist(avgDistanceKNK[gSets])
  netDist <- unlist(avgDistanceNET[gSets])
  
  expDist <- data.frame(dist = expDist, method = 'Expression')
  knkDist <- data.frame(dist = knkDist, method = 'scTenifoldKnk')
  netDist <- data.frame(dist = netDist, method = 'Network')
  allDist <- rbind(expDist, knkDist, netDist)
  allDist$method <- factor(allDist$method, levels = c('Expression', 'Network', 'scTenifoldKnk'))
  
  kTest <- kruskal.test(allDist$dist, allDist$method)
  kTest <- formatC(kTest$p.value, digits = 3)
  
  ggplot(allDist, aes(method, dist, color = method)) + geom_boxplot(outlier.color = NA) + 
    geom_jitter(alpha = 0.1, cex = 0.2) + 
    theme_bw() + 
    xlab('Method') + 
    ylab(IQR~(Mahalinobis~Distance)) + 
    #labs(title = title, subtitle = paste0('U test P-value = ', kTest)) + 
    theme(plot.title = element_text(face = 2, size = 25), legend.position = 'none', axis.text.x=element_blank()) +
    guides(color=guide_legend(title="Method:"))
}
makeBP <- function(db, title){
  tsnePositions <- S#@reductions$tsne@cell.embeddings
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  avgDistanceKNK <- lapply(db, function(gSet){
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    if(length(gSet) > 2){
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      mean(gSet)  
    }
  })
  
  tsnePositions <- eMatrix#@reductions$tsne@cell.embeddings
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  avgDistanceEXP <- lapply(db, function(gSet){
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    if(length(gSet) > 2){
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      mean(gSet)  
    }
  })
  
  tsnePositions <- aMatrix#@reductions$tsne@cell.embeddings
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  avgDistanceNET <- lapply(db, function(gSet){
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    if(length(gSet) > 2){
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      mean(gSet)  
    }
  })
  
  gSets <- intersect(names(avgDistanceEXP), names(avgDistanceKNK))
  gSets <- intersect(gSets, names(avgDistanceNET))
  
  
  expDist <- unlist(avgDistanceEXP[gSets])
  knkDist <- unlist(avgDistanceKNK[gSets])
  netDist <- unlist(avgDistanceNET[gSets])
  
  expDist <- data.frame(dist = expDist, method = 'Expression')
  knkDist <- data.frame(dist = knkDist, method = 'scTenifoldKnk')
  netDist <- data.frame(dist = netDist, method = 'Network')
  allDist <- rbind(expDist, knkDist, netDist)
  allDist$method <- factor(allDist$method, levels = c('Expression', 'Network', 'scTenifoldKnk'))
  
  kTest <- kruskal.test(allDist$dist, allDist$method)
  kTest <- formatC(kTest$p.value, digits = 3)
  
  ggplot(allDist, aes(method, dist, color = method)) + geom_boxplot(outlier.color = NA) + 
    geom_jitter(alpha = 0.1, cex = 0.2) + 
    theme_bw() + 
    xlab('Method') + 
    ylab(Avg~(Mahalinobis~Distance)) + 
    #labs(title = title, subtitle = paste0('U test P-value = ', kTest)) + 
    theme(plot.title = element_text(face = 2, size = 25), legend.position = 'none', axis.text.x=element_blank()) +
    guides(color=guide_legend(title="Method:"))
}

png('iqrkegg.png', width = 1600, height = 1000, res = 300)
A <- makeiqrBP(KEGG, 'KEGG')
B <- makeiqrECDF(KEGG, 'KEGG')
A + B + plot_layout(design = c('ABBBBB'))
dev.off()

png('kegg.png', width = 1600, height = 1000, res = 300)
A <- makeBP(KEGG, 'KEGG')
B <- makeECDF(KEGG, 'KEGG')
A + B + plot_layout(design = c('ABBBBB'))
dev.off()

png('iqrgobp.png', width = 1600, height = 1000, res = 300)
A <- makeiqrBP(GOBP, 'GO BP')
B <- makeiqrECDF(GOBP, 'GO BP')
A + B + plot_layout(design = c('ABBBBB'))
dev.off()

png('gobp.png', width = 1600, height = 1000, res = 300)
A <- makeBP(GOBP, 'GO BP')
B <- makeECDF(GOBP, 'GO BP')
A + B + plot_layout(design = c('ABBBBB'))
dev.off()

png('biop.png', width = 1600, height = 1000, res = 300)
A <- makeBP(BIOP, 'BioPlanet')
B <- makeECDF(BIOP, 'BioPlanet')
A + B + plot_layout(design = c('ABBBBB'))
dev.off()

png('iqrbiop.png', width = 1600, height = 1000, res = 300)
A <- makeiqrBP(BIOP, 'BioPlanet')
B <- makeiqrECDF(BIOP, 'BioPlanet')
A + B + plot_layout(design = c('ABBBBB'))
dev.off()

A <- FindNeighbors(S, k.param = 20)
B <- graph_from_adjacency_matrix(A$nn, weighted = TRUE)
B <- cluster_walktrap(B, steps = 50)

TSNE <- data.frame(S,B$membership)
colnames(TSNE) <- c('tSNE1', 'tSNE2', 'Cluster')
TSNE$Cluster <- as.factor(TSNE$Cluster)
TSNE$G <- rownames(S)
writeLines(TSNE$G[TSNE$Cluster %in% 4])
C2 <- TSNE[TSNE$Cluster %in% 2,]
png('knk.png', width = 1000, height = 1000, res = 300)
ggplot(TSNE, aes(tSNE1, tSNE2)) + 
  geom_point(pch = ifelse(TSNE$Cluster %in% c(2,4), 8, 16), cex = 0.5, col = ifelse(TSNE$Cluster %in% c(2,4), 'red', 'black')) + 
  labs(title = 'scTenifoldKnk') + 
  xlab('t-SNE 1') + 
  ylab('t-SNE 2') +
  #geom_encircle(data = C2, aes(tSNE1, tSNE2), lty = 2, pch = 0.5, col = 'red', lwd=2)+
  theme_bw() + 
  theme(legend.position = 'none', plot.title = element_text(face = 2, size = 15))
dev.off()


eMatrix <- eMatrix[rownames(S),]
TSNE <- data.frame(eMatrix,B$membership)
colnames(TSNE) <- c('tSNE1', 'tSNE2', 'Cluster')
TSNE$Cluster <- as.factor(TSNE$Cluster)
TSNE$G <- rownames(S)
writeLines(TSNE$G[TSNE$Cluster %in% 2])
C2 <- TSNE[TSNE$Cluster %in% 2,]
png('expression.png', width = 1000, height = 1000, res = 300)
ggplot(TSNE, aes(tSNE1, tSNE2)) + 
  geom_point(pch = ifelse(TSNE$Cluster %in% 2, 8, 16), cex = 0.5, col = ifelse(TSNE$Cluster %in% 2, 'red', 'black')) + 
  labs(title = 'Expression') + 
  xlab('t-SNE 1') + 
  ylab('t-SNE 2') +
  #geom_encircle(data = C2, aes(tSNE1, tSNE2), lty = 2, pch = 0.5, col = 'red', lwd=2)+
  theme_bw() + 
  theme(legend.position = 'none', plot.title = element_text(face = 2, size = 15))
dev.off()


aMatrix <- aMatrix[rownames(S),]
TSNE <- data.frame(aMatrix,B$membership)
colnames(TSNE) <- c('tSNE1', 'tSNE2', 'Cluster')
TSNE$Cluster <- as.factor(TSNE$Cluster)
TSNE$G <- rownames(S)
writeLines(TSNE$G[TSNE$Cluster %in% 2])
C2 <- TSNE[TSNE$Cluster %in% 2,]
png('network.png', width = 1000, height = 1000, res = 300)
ggplot(TSNE, aes(tSNE1, tSNE2)) + 
  geom_point(pch = ifelse(TSNE$Cluster %in% 2, 8, 16), cex = 0.5, col = ifelse(TSNE$Cluster %in% 2, 'red', 'black')) + 
  labs(title = 'Network') + 
  xlab('t-SNE 1') + 
  ylab('t-SNE 2') +
  #geom_encircle(data = C2, aes(tSNE1, tSNE2), lty = 2, pch = 0.5, col = 'red', lwd=2)+
  theme_bw() + 
  theme(legend.position = 'none', plot.title = element_text(face = 2, size = 15))
dev.off()
