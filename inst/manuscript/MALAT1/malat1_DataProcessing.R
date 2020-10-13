library(Matrix)
library(Seurat)
library(fgsea)
library(UpSetR)
library(harmony)
library(OrderedList)
library(igraph)
library(enrichR)
library(ggplot2)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')

mmuKEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
BIOP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
goBP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018')
WP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2019_Human')

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
DE <- FindMarkers(ALL, ident.1 = 'KO', ident.2 = 'WT', test.use = 'MAST', logfc.threshold = 0)

library(ggrepel)
DF <- data.frame(FC = DE$avg_logFC, P = -log10(DE$p_val_adj), G = rownames(DE))
DF$G[DE$p_val_adj >= 0.05] <- NA
DF$G[abs(DF$FC) < 0.29] <- NA
gCol <- densCols(DF[,1:2])
gCol[!is.na(DF$G)] <- 'red'
png('volcanoMalat1.png', width = 1200, height = 1200, res = 300)
ggplot(DF, aes(FC,P, label = G)) + 
  geom_point(col = gCol) + 
  xlim(c(-1,0.75)) + 
  geom_text_repel() +
  theme_bw() + 
  xlab(expression(log[2]~(Fold-Change))) +
  ylab(expression(-log[10]~(P-value)))
dev.off()



FC <- DE$avg_logFC #* -log10(DE$p_val_adj)
names(FC) <- toupper(rownames(DE))

#### DR ####
load('betaMALATko.RData')
write.csv(MALAT1$diffRegulation, 'drMALAT1.csv')
MALAT1$diffRegulation$gene[MALAT1$diffRegulation$p.adj < 0.05]

#### EGO PLOT #####
png('egoMalat1.png', width = 6500,height = 5000, res = 300, pointsize = 20, bg = NA)
plotKO(MALAT1, 'Malat1')
dev.off()


DR <- MALAT1$diffRegulation
Z <- DR$Z
names(Z) <- toupper(DR$gene)
writeLines(DR$gene[DR$p.adj < 0.05])


#### Comparison ####
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
eDE <- eDE[eDE$NES > 0 & eDE$padj < 0.05,]

png('upset.png', width = 750, height = 1500, res = 300)
upset(fromList(list(DE=eDE$pathway, DR=eDR$pathway)))
dev.off()

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
k <- k[k$padj < 0.05,]

png('gsea1_Malat1.png', width = 1000, height = 1000, res = 300)
NES <- round(k[k$pathway %in% 'Insulin secretion',]$NES,2)
FDR <- formatC(k[k$pathway %in% 'Insulin secretion',]$padj, format = 'e', digits = 2)
plotEnrichment(mmuKEGG$`Insulin secretion`, FC) + 
  xlab('Gene rank') + 
  ylab('Enrichment Score') + 
  theme(plot.title = element_text(face = 2)) +
  labs(title = 'Insulin secretion', subtitle = paste0('NES= ', NES, ', FDR = ', FDR))
dev.off()

png('gsea2_Malat1.png', width = 1000, height = 1000, res = 300)
NES <- round(k[k$pathway %in% 'cellular respiration (GO:0045333)',]$NES,2)
FDR <- formatC(k[k$pathway %in% 'cellular respiration (GO:0045333)',]$padj, format = 'e', digits = 2)
plotEnrichment(goBP$`cellular respiration (GO:0045333)`, FC) + 
  xlab('Gene rank') + 
  ylab('Enrichment Score') +
  theme(plot.title = element_text(face = 2)) +
  labs(title = 'Cellular Respiration GO:0045333', subtitle = paste0('NES= ', NES, ', FDR = ', FDR))
dev.off()

png('gsea3_Malat1.png', width = 1000, height = 1000, res = 300)
NES <- round(k[k$pathway %in% 'Circadian rhythm related genes WP3594',]$NES,2)
FDR <- formatC(k[k$pathway %in% 'Circadian rhythm related genes WP3594',]$padj, format = 'e', digits = 2)
plotEnrichment(WP$`Circadian rhythm related genes WP3594`, FC) + 
  xlab('Gene rank') + 
  ylab('Enrichment Score') + 
  theme(plot.title = element_text(face = 2)) +
  labs(title = 'Circadian rhythm related genes', subtitle = paste0('NES= ', NES, ', FDR = ', FDR))
dev.off()
