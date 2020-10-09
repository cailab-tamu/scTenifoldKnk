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
# SRS3161261 <- readMM('CFTR/Data/matrix_SRS3161261.mtx')
# rownames(SRS3161261) <- readLines('CFTR/Data/geneNames_SRS3161261.txt')
# colnames(SRS3161261) <- readLines('CFTR/Data/barcodes_SRS3161261.txt')
# SRS3161261 <- scQC(SRS3161261)
# SRS3161261 <- SRS3161261[rowMeans(SRS3161261 != 0) > 0.05,]
# SRS3161261 <- SRS3161261[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(SRS3161261), ignore.case = TRUE),]
# SRS3161261 <- scTenifoldKnk(SRS3161261, gKO = 'Cftr')
# save(SRS3161261, file = 'CFTR_SRS3161261.RData')
#
# SRS4245406 <- readMM('CFTR/Data/matrix_SRS4245406.mtx')
# rownames(SRS4245406) <- readLines('CFTR/Data/geneNames_SRS4245406.txt')
# colnames(SRS4245406) <- readLines('CFTR/Data/barcodes_SRS4245406.txt')
# SRS4245406 <- scQC(SRS4245406)
# SRS4245406 <- SRS4245406[rowMeans(SRS4245406 != 0) > 0.05,]
# SRS4245406 <- SRS4245406[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(SRS4245406), ignore.case = TRUE),]
# SRS4245406 <- scTenifoldKnk(SRS4245406, gKO = 'Cftr')
# save(SRS4245406, file = 'CFTR_SRS4245406.RData')

library(fgsea)
library(ggplot2)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')

load('../Results/SRS4245406.RData')

write.csv(SRS4245406$diffRegulation, file = '../Results/scTenifoldKnk_SRS4245406.csv', row.names = FALSE, quote = FALSE)
png('../Results/drSRS4245406.png', width = 2000, height = 2000, res = 300, pointsize = 20)
plotDR(SRS4245406)
dev.off()

Z <- SRS4245406$diffRegulation$Z
names(Z) <- toupper(SRS4245406$diffRegulation$gene)

CC <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2018')
BP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018')
MF <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Molecular_Function_2018')
GO <- c(CC, BP, MF)
set.seed(1)
E <- fgseaMultilevel(GO, Z)
E <- E[E$NES > 0 & E$padj < 0.05,]
E$leadingEdge <- unlist(lapply(E$leadingEdge, function(X){paste0(X, collapse = ';')}))

png('../Results/mp1_SRS4245406.png', width = 1000, height = 1000, res = 300)
gSet <- 'positive regulation of ion transmembrane transporter activity (GO:0032414)'
plotEnrichment(pathway = GO[[gSet]], stats = Z) +
  labs(title = 'Regulation of\nion transmemb.\ntransporter\n(GO:0032414)', 
       subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], format = 'e', digits = 3))) +
  xlab('Gene rank') + 
  ylab('Enrichment Score') +
  theme(plot.title = element_text(face = 2, size = 25))
dev.off()

E <- E[order(E$padj),]
write.csv('../Results/goEnrichment_CFTR.csv')

MGI <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=MGI_Mammalian_Phenotype_Level_4_2019')
E <- fgseaMultilevel(MGI, Z)
E <- E[E$padj < 0.05,]
E <- E[order(E$pval),]
E$leadingEdge <- unlist(lapply(E$leadingEdge, function(X){paste0(X, collapse = ';')}))
write.csv(E, file = '../Results/enrichmentMP.csv', row.names = FALSE, quote = FALSE)

# png('../Results/mp1_SRS4245406.png', width = 1000, height = 1000, res = 300)
# gSet <- 'MP:0004782 abnormal surfactant physiology'
# plotEnrichment(pathway = MGI[[gSet]], stats = Z) +
#   labs(title = 'Abnormal surfactant physiology\nMP:0004782', subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], format = 'e', digits = 3))) +
#   xlab('Gene rank') + ylab('Enrichment Score')
# dev.off()

png('../Results/mp2_SRS4245406.png', width = 1000, height = 1000, res = 300)
gSet <- 'MP:0004780 abnormal surfactant secretion'
plotEnrichment(pathway = MGI[[gSet]], stats = Z) +
  labs(title = 'Abnormal\nsurfactant\nsecretion\n(MP:0004780)', subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], format = 'e', digits = 3))) +
  xlab('Gene rank') + ylab('Enrichment Score') +
  theme(plot.title = element_text(face = 2, size = 25))
dev.off()

png('../Results/mp3_SRS4245406.png', width = 1000, height = 1000, res = 300)
gSet <- 'MP:0002270 abnormal pulmonary alveolus morphology'
plotEnrichment(pathway = MGI[[gSet]], stats = Z) +
  labs(title = 'Abnormal\nalveolus\nmorphology\n(MP:0002270)', subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], format = 'e', digits = 3))) +
  xlab('Gene rank') + ylab('Enrichment Score') +
  theme(plot.title = element_text(face = 2, size = 25))
dev.off()

png('../Results/ego_SRS4245406.png', width = 3000, height = 3000, res = 300, pointsize = 20)
plotKO(SRS4245406, gKO = 'Cftr', nCategories = 1)
dev.off()


