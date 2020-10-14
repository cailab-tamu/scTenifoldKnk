library(Matrix)
library(Seurat)
library(fgsea)
library(UpSetR)
library(harmony)
library(OrderedList)
library(igraph)
library(enrichR)
library(ggrepel)
library(ggplot2)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')

mmuKEGG <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
BIOP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
REACTOME <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
goBP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018')
WP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2019_Human')

#### DE ####
WT <- readMM('WT.mtx')
rownames(WT) <- readLines('genesWT.txt')
colnames(WT) <- readLines('barcodesWT.txt')
#WT <- WT[!rownames(WT) %in% 'Malat1',]

KO <- readMM('KO.mtx')
rownames(KO) <- readLines('genesKO.txt')
colnames(KO) <- readLines('barcodesKO.txt')
#KO <- KO[!rownames(KO) %in% 'Malat1',]

WT <- CreateSeuratObject(WT, project = 'WT')
KO <- CreateSeuratObject(KO, project = 'KO')

ALL <- merge(WT,KO)
ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)
ALL <- FindVariableFeatures(ALL, verbose = FALSE)
ALL <- RunPCA(ALL, verbose = FALSE)
ALL <- RunHarmony(ALL, group.by.vars = 'orig.ident')
ALL <- RunUMAP(ALL, reduction = 'harmony', dims = 1:50, verbose = FALSE)
png('umapMalat1.png', width = 1500, height = 1200, res = 300)
UMAPPlot(ALL) + theme_bw() + xlab('UMAP 1') + ylab('UMAP 2')
dev.off()



#### Volcano ####
load('betaMALATko.RData')
DE <- FindMarkers(ALL, ident.1 = 'KO', ident.2 = 'WT', test.use = 'MAST', logfc.threshold = 0)
drGenes <- MALAT1$diffRegulation$gene[MALAT1$diffRegulation$p.adj < 0.05]

library(ggrepel)
DF <- data.frame(FC = DE$avg_logFC, P = -log10(DE$p_val_adj), G = rownames(DE))
DF$G[DE$p_val_adj >= 0.05] <- NA
DF$G[abs(DF$FC) < 0.25] <- NA
#DF$G[!DF$G %in% MALAT1$diffRegulation$gene[MALAT1$diffRegulation$p.adj < 0.05]] <- NA
gCol <- densCols(DF[,1:2])
gCol[!is.na(DF$G)] <- 'red'

png('volcanoMalat1.png', width = 1200, height = 1200, res = 300)
ggplot(DF, aes(FC,P, label = G)) + 
  geom_point(col = gCol) + 
  xlim(c(-1,1)) + 
  geom_text_repel(size = 3, segment.alpha = 0.2, force = 10, segment.size = 0.1, aes(fontface = ifelse(DF$G %in% drGenes, 2, 1))) +
  theme_bw() + 
  xlab(expression(log[2]~(Fold-Change))) +
  ylab(expression(-log[10]~(P-value)))
dev.off()


#### DR ####
load('betaMALATko.RData')
write.csv(MALAT1$diffRegulation, 'drMALAT1.csv')
DR <- MALAT1$diffRegulation

# totalGenes <- length(unique(c(rownames(DE), rownames(DR))))
# deGenes <- sum(DE$p_val_adj < 0.05 & abs(DE$avg_logFC) > 0.1)
# drGenes <- sum(DR$p.adj < 0.05)
# sGenes <- length(intersect(rownames(DE)[DE$p_val_adj < 0.05 & abs(DE$avg_logFC) > 0.1], DR$gene[DR$p.adj < 0.05]))
# nrow(DR)

#### EGO PLOT #####
png('egoMalat1.png', width = 6500,height = 5000, res = 300, pointsize = 20, bg = NA)
plotKO(MALAT1, 'Malat1')
dev.off()

#### Comparison ####

Z <- DR$Z
names(Z) <- toupper(DR$gene)
writeLines(DR$gene[DR$p.adj < 0.05])
FC <- DE$avg_logFC #* -log10(DE$p_val_adj)
names(FC) <- toupper(rownames(DE))
sGenes <- intersect(names(Z), names(FC))
Z <- Z[sGenes]
FC <- abs(FC[sGenes])

plot(Z,FC)
cor(Z,FC, method = 'sp')

set.seed(1)
eDR <- fgseaMultilevel(mmuKEGG, Z)
set.seed(1)
eDE <- fgseaMultilevel(mmuKEGG, FC)

eDR <- eDR[eDR$NES > 0 & eDR$padj < 0.05,]
eDR <- eDR[order(eDR$NES, decreasing = TRUE),]
eDR$leadingEdge <- unlist(lapply(eDR$leadingEdge, function(X){paste0(X,collapse = ';')}))
write.csv(eDR, 'drEnrichment.csv')
eDE <- eDE[eDE$NES > 0 & eDE$padj < 0.05,]
eDE <- eDE[order(eDE$NES, decreasing = TRUE),]
eDE$leadingEdge <- unlist(lapply(eDE$leadingEdge, function(X){paste0(X,collapse = ';')}))
write.csv(eDE, 'deEnrichment.csv')



png('upset.png', width = 900, height = 900, res = 300, pointsize = 1.2)
upset(fromList(list(DE=eDE$pathway, DR=eDR$pathway)),mb.ratio = c(0.8, 0.2))
dev.off()


phyper(m = 85, n = 303-85, q = 8, k = 11, lower.tail = FALSE)


lZ <- names(sort(Z, decreasing = TRUE))
lFC <- names(sort(FC, decreasing = TRUE))

# png('oList.png', width = 1500, height = 1000, res = 300)
# plot(compareLists(lZ,lFC, no.reverse = TRUE, alphas = c(0.005)))
# dev.off()

allDB <- c(mmuKEGG, BIOP, WP, REACTOME, goBP)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))
#FC <- FC[DE$p_val_adj < 0.05]
FC <- FC[FC > 0]
FC <- FC[!grepl('Rpl|Rps|Rp[[:digit:]]+|Mt-',names(FC))]

set.seed(1)
kA <- fgsea(mmuKEGG, FC, nperm = 1e6)
set.seed(1)
kB <- fgsea(BIOP, FC, nperm = 1e6)
set.seed(1)
kC <- fgsea(WP, FC, nperm = 1e6)
set.seed(1)
kD <- fgsea(REACTOME, FC, nperm = 1e6)
set.seed(1)
kE <- fgsea(goBP, FC, nperm = 1e6)
k <- rbind(kA,kB,kC,kD,kE)
k1 <- k[k$padj < 0.05,]
k1$leadingEdge <- unlist(lapply(k1$leadingEdge, function(X){paste0(X, collapse = ';')}))
write.csv(k1, 'deEnrichment.csv')

set.seed(1)
kA <- fgsea(mmuKEGG, Z, nperm = 1e6)
set.seed(1)
kB <- fgsea(BIOP, Z, nperm = 1e6)
set.seed(1)
kC <- fgsea(WP, Z, nperm = 1e6)
set.seed(1)
kD <- fgsea(REACTOME, Z, nperm = 1e6)
set.seed(1)
kE <- fgsea(goBP, Z, nperm = 1e6)
k <- rbind(kA,kB,kC,kD,kE)
k2 <- k[k$padj < 0.05,]
k2$leadingEdge <- unlist(lapply(k2$leadingEdge, function(X){paste0(X, collapse = ';')}))
write.csv(k2, 'drEnrichment.csv')

k <- k2
png('gsea1_Malat1.png', width = 1000, height = 1000, res = 300)
NES <- round(k[k$pathway %in% 'Insulin secretion',]$NES,2)
FDR <- formatC(k[k$pathway %in% 'Insulin secretion',]$padj, format = 'e', digits = 2)
plotEnrichment(mmuKEGG$`Insulin secretion`, FC) + 
  xlab('Gene rank') + 
  ylab('Enrichment Score') + 
  theme(plot.title = element_text(face = 2, size = 25)) +
  labs(title = 'Insulin\nsecretion', subtitle = paste0('FDR = ', FDR))
dev.off()

# png('gsea2_Malat1.png', width = 1000, height = 1000, res = 300)
# NES <- round(k[k$pathway %in% 'cellular respiration (GO:0045333)',]$NES,2)
# FDR <- formatC(k[k$pathway %in% 'cellular respiration (GO:0045333)',]$padj, format = 'e', digits = 2)
# plotEnrichment(goBP$`cellular respiration (GO:0045333)`, FC) + 
#   xlab('Gene rank') + 
#   ylab('Enrichment Score') +
#   theme(plot.title = element_text(face = 2)) +
#   labs(title = 'Cellular Respiration GO:0045333', subtitle = paste0('NES= ', NES, ', FDR = ', FDR))
# dev.off()

png('gsea3_Malat1.png', width = 1000, height = 1000, res = 300)
NES <- round(k[k$pathway %in% 'Circadian rhythm related genes WP3594',]$NES,2)
FDR <- formatC(k[k$pathway %in% 'Circadian rhythm related genes WP3594',]$padj, format = 'e', digits = 2)
plotEnrichment(WP$`Circadian rhythm related genes WP3594`, FC) + 
  xlab('Gene rank') + 
  ylab('Enrichment Score') + 
  theme(plot.title = element_text(face = 2, size = 25)) +
  labs(title = 'Circadian rhythm\nrelated genes', subtitle = paste0('FDR = ', FDR))
dev.off()
