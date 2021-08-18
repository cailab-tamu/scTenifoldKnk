library(Matrix)
library(scTenifoldNet)
library(pbapply)
library(UpSetR)
library(ComplexHeatmap)

X <- readMM('CFTR/Data/matrix_SRS4245406.mtx')
rownames(X) <- readLines('CFTR/Data/geneNames_SRS4245406.txt')
colnames(X) <- readLines('CFTR/Data/barcodes_SRS4245406.txt')

X <- X[rowMeans(X != 0) > 0.05,]
X <- as.matrix(X)

spCor <- pbsapply(rownames(X), function(G){
  cor(X['Cftr',], X[G,], method = 'sp')
})

spCor <- sort(abs(spCor), decreasing = TRUE)

write.csv(spCor, file = 'reviewer3_comment2_corr.csv')


highlyCorrelated <- names(spCor)[2:11]

load('CFTR/Results/SRS4245406.RData')

WT <- as.matrix(SRS4245406$WT)
KO <- WT

drGenes <- pblapply(seq_along(highlyCorrelated), function(i){
  KO[highlyCorrelated[i],] <- 0
  MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2)
  DR <- scTenifoldKnk:::dRegulation(MA, highlyCorrelated[i])
  write.csv(DR, paste0('reviewer3_comment2_ko_',highlyCorrelated[i],'.csv'))
  DR$gene[DR$p.adj < 0.05]  
})

names(drGenes) <- highlyCorrelated
upset(fromList(drGenes), nsets = 10, nintersects = 100)

allDR <- table(unlist(drGenes))
writeLines(names(allDR[allDR >= 10]), sep = ', ')

allZ <- lapply(highlyCorrelated, function(i){
  DR <- read.csv(paste0('reviewer3_comment2_ko_',i,'.csv'))
  Z <- DR$Z
  names(Z) <- DR$gene
  return(Z)
})

allZ <- sapply(allZ, function(X){X[unique(unlist(lapply(allZ, names)))]})
colnames(allZ) <- highlyCorrelated
mean(cor(allZ, method = 'sp'))
sd(cor(allZ, method = 'sp'))

alwaysDR <- names(allDR[allDR >= 10])
at <- which(rownames(allZ) %in% alwaysDR)
label <- rownames(allZ)[at]
png('reviewer3_comment2.png', width = 1000, height = 2000, res = 300)
Heatmap(allZ, show_row_names = FALSE, show_row_dend = FALSE, name = 'Z') + 
  rowAnnotation(link = anno_mark(at = at, labels = label))
dev.off()
