# setwd('/data/dcosorioh/manuscript/')
#
# library(Matrix)
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
# source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')
#
# WT <- readMM('TREM2/Data/matrix_WT.mtx')
# rownames(WT) <- readLines('TREM2/Data/geneNames_WT.txt')
# colnames(WT) <- readLines('TREM2/Data/barcodes_WT.txt')
# WT <- scQC(WT)
# WT <- WT[rowMeans(WT != 0) > 0.05,]
# WT <- WT[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(WT), ignore.case = TRUE),]
# WT <- scTenifoldKnk(WT, gKO = 'Trem2')
# save(WT, file = 'TREM2_GSE130626.RData')

library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa2mmu_SYMBOL.R')

load('../Results/GSE130626.RData')
strictDirection <- function(X){
  lT <- X[lower.tri(X)]
  uT <- t(X)[lower.tri(X)]
  lT[abs(lT) < abs(uT)] <- 0
  uT[abs(uT) < abs(lT)] <- 0
  X[lower.tri(X)] <- lT
  X[upper.tri(X)] <- uT
  X <- Matrix::Matrix(X)
  return(X)
}

library(Matrix)
X <- GSE130626$WT
X <- strictDirection(X)

# sum((GSE130626$WT - X) != 0)
# Y <- X
# Y['Trem2',] <- 0
# MA <- scTenifoldNet::manifoldAlignment(X,Y,2)
# DR <- scTenifoldKnk:::dRegulation(MA[!grepl('_RPL|_RPS|_MT-|_RP[[:digit:]]',rownames(MA), ignore.case = TRUE),], 1)
# write.csv(DR, '../dirResults/dirTREM2.csv')
#
# Y <- X
# Y[,'Trem2'] <- 0
# MA <- scTenifoldNet::manifoldAlignment(X,Y,2)
# DR <- scTenifoldKnk:::dRegulation(MA[!grepl('_RPL|_RPS|_MT-|_RP[[:digit:]]',rownames(MA), ignore.case = TRUE),], 1)
# write.csv(DR, '../dirResults/dir_colTREM2.csv')
#
# Y <- X
# Y[,'Trem2'] <- 0
# MA <- scTenifoldNet::manifoldAlignment(t(X),t(Y),2)
# DR <- scTenifoldKnk:::dRegulation(MA[!grepl('_RPL|_RPS|_MT-|_RP[[:digit:]]',rownames(MA), ignore.case = TRUE),], 1)
# write.csv(DR, '../dirResults/dir_coltransposeTREM2.csv')


dirDR <- read.csv('../dirResults/dir_coltransposeTREM2.csv', row.names = 1)
oriDR <- GSE130626$diffRegulation

Z1 <- dirDR$Z
names(Z1) <- toupper(dirDR$gene)

Z2 <- oriDR$Z
names(Z2) <- toupper(oriDR$gene)
iGenes <- intersect(names(Z1), names(Z2))
cor(Z2[iGenes],Z1[iGenes], method = 'sp')
plot(Z2[iGenes],Z1[iGenes])

GSE130626$diffRegulation <- dirDR
GSE130626$WT <- t(X)
plotKO(GSE130626, gKO = 'Trem2')

dGenes <- GSE130626$diffRegulation$gene[GSE130626$diffRegulation$p.adj < 0.05]

png('../dirResults/colGSE130626.png', width = 2000, height = 2000, res = 300)
plotDR(GSE130626)
dev.off()

png('../dirResults/ego_colGSE130626.png', width = 3000, height = 3000, res = 300, bg = NA)
plotKO(GSE130626, gKO = 'Trem2', nCategories = 4)
dev.off()


zGSE130626 <- GSE130626$diffRegulation$Z
names(zGSE130626) <- toupper(GSE130626$diffRegulation$gene)

MGI <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
set.seed(1)
E <- fgseaMultilevel(MGI, zGSE130626)

png('../dirResults/gsea1_colGSE130626.png', width = 1000, height = 1000, res = 300)
gSet <- 'Cholesterol metabolism'
plotEnrichment(MGI[[gSet]], zGSE130626) +
  labs(
    title = 'Cholesterol\nmetabolism',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

png('../dirResults/gsea2_colGSE130626.png', width = 1000, height = 1000, res = 300)
gSet <- 'Alzheimer disease'
plotEnrichment(MGI[[gSet]], zGSE130626) +
  labs(
    title = 'Alzheimer\ndisease',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

png('../dirResults/gsea3_colGSE130626.png', width = 1000, height = 1000, res = 300)
gSet <- 'Lysosome'
plotEnrichment(MGI[[gSet]], zGSE130626) +
  labs(
    title = gSet,
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score')  + theme(plot.title = element_text(face = 2, size = 25))
dev.off()
