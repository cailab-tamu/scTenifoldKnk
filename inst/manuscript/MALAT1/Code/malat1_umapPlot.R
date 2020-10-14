library(Matrix)
library(Seurat)
library(harmony)
library(fgsea)
library(ggplot2)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')

WT <- Read10X_h5('WT.h5')
WT <- scQC(WT, mtThreshold = 0.05)
WT <- CreateSeuratObject(WT, project = 'WT')
KO <- Read10X_h5('KO.h5')
KO <- scQC(KO, mtThreshold = 0.05)
KO <- CreateSeuratObject(KO, project = 'KO')
ALL <- merge(WT,KO)
ALL <- NormalizeData(ALL)
ALL <- FindVariableFeatures(ALL)
ALL <- ScaleData(ALL)
ALL <- RunPCA(ALL)
ALL <- RunHarmony(ALL, group.by.vars = 'orig.ident')
ALL <- RunUMAP(ALL, reduction = 'harmony', dims = 1:20)
UMAPPlot(ALL)
ALL <- RunTSNE(ALL, reduction = 'harmony', dims = 1:20)
TSNEPlot(ALL)
ALL <- FindNeighbors(ALL, reduction = 'tsne', dims = 1:2)
ALL <- FindClusters(ALL, resolution = 0.01)
UMAPPlot(ALL)
ctDE <- FindAllMarkers(ALL, min.pct = 0.05)

CT <- gmtPathways('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/markerGenes/mmuPanglaoDB.gmt')

ctList <- sapply(levels(Idents(ALL)), function(cID){
  FC <- ctDE$avg_logFC[ctDE$cluster %in% cID]
  names(FC) <-  ctDE$gene[ctDE$cluster %in% cID]
  set.seed(1)
  E <- fgseaMultilevel(CT, FC)
  E <- E[order(1/abs(E$NES), E$padj),]
  E <- E[E$NES > 0 & E$padj < 0.05,]
  E$pathway[1]
})
levels(Idents(ALL)) <- ctList
png('malat1_CellTypes.png', width = 1000, height = 500, res = 300, pointsize = 5)
TSNEPlot(ALL) + theme_bw() + xlab('t-SNE 1') + ylab('t-SNE 2')
dev.off()
UMAPPlot(ALL)

png('malat1_Class.png', width = 800, height = 500, res = 300, pointsize = 5)
Idents(ALL) <- ALL$orig.ident
TSNEPlot(ALL) + theme_bw() + xlab('t-SNE 1') + ylab('t-SNE 2')
dev.off()
UMAPPlot(ALL)
