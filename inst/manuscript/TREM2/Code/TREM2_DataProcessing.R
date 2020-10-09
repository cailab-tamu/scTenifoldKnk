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
dGenes <- GSE130626$diffRegulation$gene[GSE130626$diffRegulation$p.adj < 0.05]

png('../Results/dr2_GSE130626.png', width = 2000, height = 2000, res = 300)
plotDR(GSE130626)
dev.off()

png('../Results/ego2_GSE130626.png', width = 3000, height = 3000, res = 300, bg = NA)
X <- GSE130626
gKO <- 'Trem2'
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
E <- E[E$Term %in% c('Oxidative phosphorylation','Alzheimer disease','Cholesterol metabolism','Lysosome','Neutrophil mediated immunity'),]
E <- E[c(1,2,5,6,9),]
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
vColor <- rgb(195/255, 199/255, 198/255 ,0.3)
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


zGSE130626 <- GSE130626$diffRegulation$Z
names(zGSE130626) <- toupper(GSE130626$diffRegulation$gene)

MGI <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
set.seed(1)
E <- fgseaMultilevel(MGI, zGSE130626)

png('../Results/gsea1_GSE130626.png', width = 1000, height = 1000, res = 300)
gSet <- 'Cholesterol metabolism'
plotEnrichment(MGI[[gSet]], zGSE130626) +
  labs(
    title = 'Cholesterol\nmetabolism',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

png('../Results/gsea2_GSE130626.png', width = 1000, height = 1000, res = 300)
gSet <- 'Alzheimer disease'
plotEnrichment(MGI[[gSet]], zGSE130626) +
  labs(
    title = 'Alzheimer\ndisease',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

png('../Results/gsea3_GSE130626.png', width = 1000, height = 1000, res = 300)
gSet <- 'Lysosome'
plotEnrichment(MGI[[gSet]], zGSE130626) +
  labs(
    title = gSet,
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score')  + theme(plot.title = element_text(face = 2, size = 25))
dev.off()
