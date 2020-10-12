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
# WT <- readMM('MECP2/Data/matrix_SRS3059998.mtx')
# rownames(WT) <- readLines('MECP2/Data/geneNames_SRS3059998.txt')
# colnames(WT) <- readLines('MECP2/Data/barcodes_SRS3059998.txt')
# WT <- scQC(WT)
# WT <- WT[rowMeans(WT != 0) > 0.05,]
# WT <- WT[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(WT), ignore.case = TRUE),]
# WT <- scTenifoldKnk(WT, gKO = 'Mecp2')
# save(WT, file = 'MECP2_SRS3059998.RData')
#
# WT <- readMM('MECP2/Data/matrix_SRS3059999.mtx')
# rownames(WT) <- readLines('MECP2/Data/geneNames_SRS3059999.txt')
# colnames(WT) <- readLines('MECP2/Data/barcodes_SRS3059999.txt')
# WT <- scQC(WT)
# WT <- WT[rowMeans(WT != 0) > 0.05,]
# WT <- WT[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(WT), ignore.case = TRUE),]
# WT <- scTenifoldKnk(WT, gKO = 'Mecp2')
# save(WT, file = 'MECP2_SRS3059999.RData')
#
library(UpSetR)
library(pbapply)
library(ggplot2)
library(enrichR)
library(igraph)
library(fgsea)
library(statsExpressions)
library(patchwork)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')


load('../Results/SRS3059998.RData')
load('../Results/SRS3059999.RData')

write.csv(SRS3059998$diffRegulation, '../Results/drSRS3059998.csv')
write.csv(SRS3059999$diffRegulation, '../Results/drSRS3059999.csv')

gList1 <- SRS3059998$diffRegulation$gene[SRS3059998$diffRegulation$p.adj < 0.05]
gList2 <- SRS3059999$diffRegulation$gene[SRS3059999$diffRegulation$p.adj < 0.05]

sharedGene <- intersect(gList1, gList2)

DB <- listEnrichrDbs()
TF <- enrichr(sharedGene, databases = 'ENCODE_TF_ChIP-seq_2015')
TF <- do.call(rbind.data.frame, TF)
TF <- TF[order(TF$Adjusted.P.value),]
write.csv(TF, file = '../Results/tf_MECP2.csv', row.names = FALSE, quote = FALSE)
TF <- TF[1:10,]
dF <- data.frame(TF=TF$Term, FDR = -log10(TF$Adjusted.P.value), OVLP = TF$Overlap)
dF$TF <- factor(dF$TF, levels = rev(dF$TF))
colorList <- hcl.colors(3,'Zissou 1')
colorList <- rev(c(rep(colorList[3], 8), colorList[2:1]))

png(filename = '../Results/tf_MECP2.png',width = 750, height = 500, res = 300)
ggplot(dF, mapping = aes(x = TF, y = FDR, label= OVLP)) +
  geom_bar(stat = 'identity', fill = colorList) +
  geom_text(color="white", size=2, nudge_y = -1.4) +
  coord_flip() + theme_minimal() + xlab('ENCODE TF') + ylab(expression(-log[10] * '(FDR)'))
dev.off()

MET <- enrichr(sharedGene, c('Reactome_2016','BioPlanet_2019','KEGG_2019_Mouse', 'GO_Biological_Process_2018', 'GO_Molecular_Function_2018', 'GO_Cellular_Component_2018'))
MET <- do.call(rbind.data.frame, MET)
MET <- MET[order(MET$Adjusted.P.value),]
MET <- MET[MET$Adjusted.P.value < 0.05,]
write.csv(MET, '../Results/mp_MECP2.csv',quote = FALSE)

MET <- MET[c(1,2,3,10,15,18,19,23,24,27),]
MET$Term <- unlist(lapply(strsplit(MET$Term, ' Homo| \\('), function(X){X[1]}))
MET$Term <- unlist(lapply(strsplit(MET$Term,''), function(X){X[1] <- paste0(toupper(X[1]), paste0(X[-1],collapse = ''), collapse = '')}))
dF <- data.frame(TF=MET$Term, FDR = -log10(MET$Adjusted.P.value), OVLP = MET$Overlap)
dF$TF <- factor(dF$TF, levels = rev(dF$TF))
png(filename = '../Results/enrichment_MECP2.png',width = 1250, height = 500, res = 300)
ggplot(dF, mapping = aes(x = TF, y = FDR, label= OVLP)) +
  geom_bar(stat = 'identity', fill = hcl.colors(3,'Zissou 1')[1]) +
  geom_text(color="white", size=2, nudge_y = -1.7) +
  coord_flip() + theme_minimal() + xlab('Gene Set') + ylab(expression(-log[10] * '(FDR)'))
dev.off()

png('../Results/upset_MECP2.png', width = 1300, height = 1500, res = 300)
upset(fromList(list(SRS3059998 = gList1, SRS3059999 = gList2)), text.scale = 1.5)
dev.off()

gList1 <- SRS3059998$diffRegulation$gene
gList2 <- SRS3059999$diffRegulation$gene


overlapR <- pbsapply(seq_len(10), function(Z){
  nRanks <- min(c(length(gList1), length(gList2)))
  rList1 <- sample(gList1)
  rList2 <- sample(gList2)
  rOverlap <- sapply(seq_len(nRanks), function(X){
    length(intersect(rList1[seq_len(X)], rList2[seq_len(X)]))
  })
  return(rOverlap)
})

nRanks <- min(c(length(gList1), length(gList2)))
rOverlap <- sapply(seq_len(nRanks), function(X){
  length(intersect(gList1[seq_len(X)], gList2[seq_len(X)]))
})

listComparison <- data.frame(
  rank = seq_len(nRanks),
  real = rOverlap,
  expLB = apply(overlapR,1,mean) - 5* apply(overlapR,1,sd),
  exp = apply(overlapR,1,mean),
  expUB = apply(overlapR,1,mean) + 5* apply(overlapR,1,sd)
)
png('../Results/cl_MECP2.png', width = 1800, height = 1000, res = 300, pointsize = 20)
ggplot(data = listComparison, mapping = aes(x = rank, y = real)) +
  geom_line(mapping = aes(colour='SRS3059998\nSRS3059999')) +
  geom_ribbon(mapping = aes(ymin=expLB, ymax=expUB), alpha=0.3, fill = "#E69F00") +
  geom_line(mapping = aes(x = rank, y = exp, colour='Expected by\nrandom')) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() + xlab('Rank') + ylab('Size of Overlap') + ylim(c(0,7776)) + xlim(c(0,7776)) +
  scale_fill_manual(name="Overlap:",aesthetics = 'colour', values=c("#E69F00","#999999")) +
  theme(legend.position="top", axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18))
dev.off()


X <- SRS3059998
gKO <- 'Mecp2'
q <- 1
library(enrichR)
library(igraph)
gList <- unique(c(gKO, X$diffRegulation$gene[X$diffRegulation$distance > 1e-10 & X$diffRegulation$p.adj < 0.05]))
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
uniqueColor <- hcl.colors(5, palette = 'Zissou 1', alpha = 0.5)[1]
sharedColor <- hcl.colors(5, palette = 'Zissou 1', alpha = 0.5)[5]
vColor <- ifelse(names(V(netPlot)) %in% sharedGene, sharedColor, uniqueColor)
png('../Results/ego_SRS3059998.png', width = 6000, height = 6000, pointsize = 15, res=300, bg = NA)
par(mar=c(0,0,0,0))
plot(netPlot,
         layout = layPlot,
         edge.arrow.size=.2,
         vertex.label.color="black",
         vertex.size = 10+dPlot,
         vertex.label.family="Times",
         vertex.label.font=ifelse(names(V(netPlot)) %in% sharedGene,2,1),
         edge.color = ifelse(E(netPlot)$W > 0, 'red', 'blue'),
         edge.curved = ifelse(W == 0.2, 0, 0.1),
         vertex.color = vColor,
         vertex.frame.color = NA)
text(0,1.15, "SRS3059998", font = 2, cex = 2)
dev.off()


X <- SRS3059999
gKO <- 'Mecp2'
q <- 1
gList <- unique(c(gKO, X$diffRegulation$gene[X$diffRegulation$distance > 1e-10 & X$diffRegulation$p.adj < 0.05]))
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
uniqueColor <- hcl.colors(5, palette = 'Zissou 1', alpha = 0.5)[1]
sharedColor <- hcl.colors(5, palette = 'Zissou 1', alpha = 0.5)[5]
vColor <- ifelse(names(V(netPlot)) %in% sharedGene, sharedColor, uniqueColor)
png('../Results/ego_SRS3059999.png', width = 6000, height = 6000, pointsize = 15, res=300, bg = NA)
par(mar=c(0,0,0,0))
plot(netPlot,
     layout = layPlot,
     edge.arrow.size=.2,
     vertex.label.color="black",
     vertex.size = 10+dPlot,
     vertex.label.family="Times",
     vertex.label.font=ifelse(names(V(netPlot)) %in% sharedGene,2,1),
     edge.color = ifelse(E(netPlot)$W > 0, 'red', 'blue'),
     edge.curved = ifelse(W == 0.2, 0, 0.1),
     vertex.color = vColor,
     vertex.frame.color = NA)
text(0,1.15, "SRS3059999", font = 2, cex = 2)
dev.off()

TF <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=ChEA_2016')

zSRS3059998 <- SRS3059998$diffRegulation$Z
names(zSRS3059998) <- toupper(SRS3059998$diffRegulation$gene)
zSRS3059999 <- SRS3059999$diffRegulation$Z
names(zSRS3059999) <- toupper(SRS3059999$diffRegulation$gene)

sharedGene <- intersect(gList1,gList2)
cor(zSRS3059998[toupper(sharedGene)], zSRS3059999[toupper(sharedGene)], method = 'sp')
zScores <- data.frame(SRS3059998 = zSRS3059998[toupper(sharedGene)], SRS3059999 = zSRS3059999[toupper(sharedGene)])

png('../Results/zScores_MECP2.png', width = 1800, height = 1500, res = 300, pointsize = 20)
ggplot(data = zScores, mapping = aes(SRS3059998, SRS3059999)) +
  geom_point(color = densCols(zScores, colramp = hcl.colors)) +
  theme_bw() +
  labs(
    subtitle = expr_corr_test(data = zScores, x = SRS3059998, y = SRS3059999, type = 'spearman')
    ) +
  geom_density_2d() +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  xlab('Z-score (SRS3059998)') +
  ylab('Z-score (SRS3059999)') +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18))
dev.off()

WP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2019_Human')
set.seed(1)
E <- fgseaMultilevel(WP,zSRS3059998)
A <- plotEnrichment(WP$`Brain-Derived Neurotrophic Factor (BDNF) signaling pathway WP2380`, stats = zSRS3059998) +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'SRS3059998\nBDNF signaling\npathway',
       subtitle = paste0('FDR = ', formatC(E[grepl('\\(BDNF\\)',E$pathway),]$padj, digits = 2, format = 'e'))) +
  theme(plot.title = element_text(face = 2, size = 25))
set.seed(1)
E <- fgseaMultilevel(WP,zSRS3059999)
B <- plotEnrichment(WP$`Brain-Derived Neurotrophic Factor (BDNF) signaling pathway WP2380`, stats = zSRS3059999) +
  xlab('Gene rank') +
  ylab('Enrichment Score') +
  labs(title = 'SRS3059999\nBDNF signaling\npathway',
       subtitle = paste0('FDR = ', formatC(E[grepl('\\(BDNF\\)',E$pathway),]$padj, digits = 2, format = 'e'))) +
  theme(plot.title = element_text(face = 2, size = 25))
png('../Results/bndf_SRS3059998.png', width = 1000, height = 750, res = 300)
A
dev.off()
png('../Results/bndf_SRS3059999.png', width = 1000, height = 750, res = 300)
B
dev.off()
