# library(pbapply)
# load('TREM2/Results/GSE130626.RData')
# WT <- GSE130626$WT
# sum(WT['Trem2', ] != 0)
# set.seed(1)
# allDR <- pblapply(1:100, function(X){
#   sGenes <- unique(c('Trem2', sample(rownames(WT), 1000)))
#   tWT <- WT[sGenes, sGenes]
#   tKO <- tWT
#   tKO['Trem2',] <- 0
#   
#   DR <- scTenifoldKnk:::dRegulation(scTenifoldNet::manifoldAlignment(tWT, tKO, d = 2),1)
#   return(DR)
# })
# drGenes <- lapply(allDR, function(X){X$gene[X$p.adj < 0.05]})
# allGenes <- lapply(allDR, function(X){X$gene})
# allGenes <- table(unlist(allGenes))
# drGenes <- table(unlist(drGenes))
# sGenes <- intersect(names(allGenes), names(drGenes))
# O <- data.frame(n.tested=(drGenes[sGenes]/allGenes[sGenes]), n = allGenes[sGenes])
# O <- O[order(O$n.tested.Freq, decreasing = TRUE),]
# writeLines(as.vector(O$n.tested.Var1)[O$n.tested.Freq > 0.8])

setwd('/data/dcosorioh/TREM2/')
library(Matrix)
library(scTenifoldKnk)

# TREM2
CM <- read.csv('Data/GSE130626_umi_counts.csv.gz')
MD <- read.csv('Data/GSE130626_cell_info.csv.gz', stringsAsFactors = FALSE)
GD <- read.csv('Data/GSE130626_gene_info.csv.gz', stringsAsFactors = FALSE)

CM <- CM[!is.na(GD$symbol),]
GD <- GD[!is.na(GD$symbol),]

rownames(CM) <- make.unique(GD$symbol)

MD <- MD[((MD$treatment %in% 'cuprizone') & (MD$trem2_genotype %in% c('WT'))),]
CM <- CM[,MD$cell_id]
CM <- Matrix(as.matrix(CM))

# set.seed(1)
# tCM <- CM[,sample(colnames(CM), size = 0.75 * ncol(CM))]
# tCM <- tCM[rowMeans(tCM != 0) > 0.1,]
# O <- scTenifoldKnk::scTenifoldKnk(tCM, gKO = 'Trem2')
# save(O, file = 'Trem2r1.RData')

# set.seed(2)
# tCM <- CM[,sample(colnames(CM), size = 0.75 * ncol(CM))]
# tCM <- tCM[rowMeans(tCM != 0) > 0.1,]
# O <- scTenifoldKnk::scTenifoldKnk(tCM, gKO = 'Trem2')
# save(O, file = 'Trem2r2.RData')

# set.seed(3)
# tCM <- CM[,sample(colnames(CM), size = 0.75 * ncol(CM))]
# tCM <- tCM[rowMeans(tCM != 0) > 0.1,]
# O <- scTenifoldKnk::scTenifoldKnk(tCM, gKO = 'Trem2')
# save(O, file = 'Trem2r3.RData')

# set.seed(4)
# tCM <- CM[,sample(colnames(CM), size = 0.75 * ncol(CM))]
# tCM <- tCM[rowMeans(tCM != 0) > 0.1,]
# O <- scTenifoldKnk::scTenifoldKnk(tCM, gKO = 'Trem2')
# save(O, file = 'Trem2r4.RData')

# set.seed(5)
# tCM <- CM[,sample(colnames(CM), size = 0.75 * ncol(CM))]
# tCM <- tCM[rowMeans(tCM != 0) > 0.1,]
# O <- scTenifoldKnk::scTenifoldKnk(tCM, gKO = 'Trem2')
# save(O, file = 'Trem2r5.RData')

# set.seed(6)
# tCM <- CM[,sample(colnames(CM), size = 0.75 * ncol(CM))]
# tCM <- tCM[rowMeans(tCM != 0) > 0.1,]
# O <- scTenifoldKnk::scTenifoldKnk(tCM, gKO = 'Trem2')
# save(O, file = 'Trem2r6.RData')

# set.seed(7)
# tCM <- CM[,sample(colnames(CM), size = 0.75 * ncol(CM))]
# tCM <- tCM[rowMeans(tCM != 0) > 0.1,]
# O <- scTenifoldKnk::scTenifoldKnk(tCM, gKO = 'Trem2')
# save(O, file = 'Trem2r7.RData')

# set.seed(8)
# tCM <- CM[,sample(colnames(CM), size = 0.75 * ncol(CM))]
# tCM <- tCM[rowMeans(tCM != 0) > 0.1,]
# O <- scTenifoldKnk::scTenifoldKnk(tCM, gKO = 'Trem2')
# save(O, file = 'Trem2r8.RData')

# set.seed(9)
# tCM <- CM[,sample(colnames(CM), size = 0.75 * ncol(CM))]
# tCM <- tCM[rowMeans(tCM != 0) > 0.1,]
# O <- scTenifoldKnk::scTenifoldKnk(tCM, gKO = 'Trem2')
# save(O, file = 'Trem2r9.RData')

# set.seed(10)
# tCM <- CM[,sample(colnames(CM), size = 0.75 * ncol(CM))]
# tCM <- tCM[rowMeans(tCM != 0) > 0.1,]
# O <- scTenifoldKnk::scTenifoldKnk(tCM, gKO = 'Trem2')
# save(O, file = 'Trem2r10.RData')