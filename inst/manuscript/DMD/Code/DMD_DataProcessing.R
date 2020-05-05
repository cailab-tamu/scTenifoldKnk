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
# GSM4116571 <- readMM('DMD/Data/GSM4116571_Qc_matrix.mtx.gz')
# rownames(GSM4116571) <- read.csv('DMD/Data/GSM4116571_Qc_features.tsv.gz', sep = '\t', header = FALSE, stringsAsFactors = FALSE)[,2]
# colnames(GSM4116571) <- readLines('DMD/Data/GSM4116571_Qc_barcodes.tsv.gz')
# GSM4116571 <- scQC(GSM4116571)
# GSM4116571 <- GSM4116571[rowMeans(GSM4116571 != 0) > 0.05,]
# GSM4116571 <- GSM4116571[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(GSM4116571), ignore.case = TRUE),]
# GSM4116571 <- scTenifoldKnk(GSM4116571, gKO = 'Dmd')
# save(GSM4116571, file = 'DMD_GSM4116571.RData')


library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa2mmu_SYMBOL.R')

load('../Results/GSM4116571.RData')
write.csv(GSM4116571$diffRegulation, file = '../Results/scTenifoldKnk_SRS4245406.csv', row.names = FALSE, quote = FALSE)

png('../Results/drGSM4116571.png', width = 2000, height = 2000, res = 300, pointsize = 20)
plotDR(GSM4116571, labelGenes = hsa2mmu_SYMBOL(readLines('../Results/enrichedGenes.txt')))
dev.off()

Z <- GSM4116571$diffRegulation$Z
names(Z) <- toupper(GSM4116571$diffRegulation$gene)

MGI <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=MGI_Mammalian_Phenotype_Level_4_2019')
E <- fgseaMultilevel(MGI, Z)
E <- E[E$padj < 0.05 & E$NES > 0,]
E <- E[order(1/E$NES,E$pval),]
E$leadingEdge <- unlist(lapply(E$leadingEdge, function(X){paste0(X, collapse = ';')}))
write.csv(E, file = '../Results/enrichmentMP.csv', row.names = FALSE, quote = FALSE)


png('../Results/mp1_GSM4116571.png', width = 1000, height = 1000, res = 300)
gSet <- 'MP:0008438 abnormal cutaneous collagen fibril morphology'
plotEnrichment(pathway = MGI[[gSet]], stats = Z) +
  labs(title = 'Abnormal collagen fibril\nmorphology MP:0008438', subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], format = 'e', digits = 3))) +
  xlab('Gene rank') + ylab('Enrichment Score')
dev.off()
png('../Results/mp2_GSM4116571.png', width = 1000, height = 1000, res = 300)
gSet <- 'MP:0000759 abnormal skeletal muscle morphology'
plotEnrichment(pathway = MGI[[gSet]], stats = Z) +
  labs(title = 'Abnormal skeletal muscle \nmorphology MP:0000759', subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], format = 'e', digits = 3))) +
  xlab('Gene rank') + ylab('Enrichment Score')
dev.off()
png('../Results/mp3_GSM4116571.png', width = 1000, height = 1000, res = 300)
gSet <- 'MP:0003084 abnormal skeletal muscle fiber morphology'
plotEnrichment(pathway = MGI[[gSet]], stats = Z) +
  labs(title = 'Abnormal skeletal muscle\nfiber morphology MP:0003084', subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], format = 'e', digits = 3))) +
  xlab('Gene rank') + ylab('Enrichment Score')
dev.off()

png('../Results/ego_GSM4116571.png', width = 3000, height = 3000, res = 300, bg = NA)
X <- GSM4116571
gKO <- 'Dmd'
q <- 0.995
gList <- unique(c(gKO, X$diffRegulation$gene[X$diffRegulation$p.adj < 0.05]))
sCluster <- as.matrix(X$WT[gList,gList])
koInfo <- sCluster[gKO,]
sCluster[abs(sCluster) <= quantile(abs(sCluster), q)] <- 0
sCluster[gKO,] <- koInfo
diag(sCluster) <- 0
sCluster <-  reshape2::melt(as.matrix(sCluster))
colnames(sCluster) <- c('from', 'to', 'W')
sCluster <- sCluster[sCluster$W != 0,]
netPlot <- graph_from_data_frame(sCluster, directed = TRUE)
dPlot <- centr_degree(netPlot)$res
W <- rep(1,nrow(sCluster))
sG   <- (names(V(netPlot))[dPlot > 1])[-1]
W[sCluster$from %in% sG] <- 0.2
W[sCluster$to %in% sG] <- 0.2
W[sCluster$from %in% gKO] <- 1
W[sCluster$from %in% gKO & sCluster$to %in% sG] <- 0.8
set.seed(1)
layPlot <- layout_with_fr(netPlot, weights = W)
dPlot <- (dPlot/max(dPlot))*20
E <- enrichr(gList, c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))
E <- do.call(rbind.data.frame, E)
E <- E[E$Adjusted.P.value < 0.05,]
E <- E[order(E$Adjusted.P.value),]
E$Term <- unlist(lapply(strsplit(E$Term,''), function(X){
  X[1] <- toupper(X[1])
  X <- paste0(X,collapse = '')
  X <- gsub('\\([[:print:]]+\\)|Homo[[:print:]]+|WP[[:digit:]]+','',X)
  X <- gsub("'s",'',X)
  X <- unlist(strsplit(X,','))[1]
  X <- gsub('[[:blank:]]$','',X)
  return(X)
}))
E <- E[E$Term %in% c("Actin filament", "Contractile actin filament bundle", "Actomyosin", "Focal adhesion", "Assembly of collagen fibrils and other multimeric structures", "Extracellular matrix organization"),]
E$Term <- gsub('Assembly of collagen fibrils and other multimeric structures','Assembly of collagen fibrils',E$Term)
E <- E[c(1,2,3,4,5,6,8),]
tPlot <- strsplit(E$Genes, ';')
pPlot <- matrix(0,nrow = length(V(netPlot)), ncol = nrow(E))
rownames(pPlot) <- toupper(names(V(netPlot)))
for(i in seq_along(tPlot)){
  pPlot[unlist(tPlot[i]),i] <- 1
}
pPlot <- lapply(seq_len(nrow(pPlot)), function(X){as.vector(pPlot[X,])})
names(pPlot) <- names(V(netPlot))
tPlot <- unique(unlist(tPlot))
eGenes <- toupper(names(V(netPlot))) %in% tPlot
vColor <- rgb(0,188/255,1,0.3)
pieColors <- list(hcl.colors(nrow(E), palette = 'Zissou 1', alpha = 0.7))
par(mar=c(4,0,0,0), xpd = TRUE)
suppressWarnings(plot(netPlot,
                      layout = layPlot,
                      edge.arrow.size=.2,
                      vertex.label.color="black",
                      vertex.shape = ifelse(eGenes,'pie','circle'),
                      vertex.pie = pPlot,
                      vertex.size = 10+dPlot,
                      vertex.pie.color=pieColors,
                      vertex.label.family="Times",
                      vertex.label.font=ifelse(eGenes,2,1),
                      edge.color = ifelse(E(netPlot)$W > 0, 'red', 'blue'),
                      edge.curved = ifelse(W == 0.2, 0, 0.1),
                      vertex.color = vColor,
                      vertex.frame.color = NA))
sigLevel <- formatC(E$Adjusted.P.value, digits = 2, format = 'g', width = 0, drop0trailing = TRUE)
gSetNames <- lengths(strsplit(E$Genes, ';'))
gSetNames <- paste0('(', gSetNames,') ', E$Term, ' FDR = ', sigLevel)
legend(x = -1.05, y = -1.05, legend = gSetNames, bty = 'n', ncol = 2, cex = 1, col = unlist(pieColors), pch = 16)
dev.off()
