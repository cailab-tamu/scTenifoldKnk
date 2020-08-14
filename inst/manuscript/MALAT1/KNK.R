library(Matrix)
library(Seurat)
library(scTenifoldKnk)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')

MALAT1 <- Read10X_h5('WT.h5')
MALAT1 <- scQC(MALAT1, mtThreshold = 0.05)
MALAT1 <- CreateSeuratObject(MALAT1)
MALAT1 <- NormalizeData(MALAT1)
MALAT1 <- FindVariableFeatures(MALAT1)
MALAT1 <- ScaleData(MALAT1)
MALAT1 <- RunPCA(MALAT1, verbose = FALSE)
MALAT1 <- RunUMAP(MALAT1, dims = 1:20)
MALAT1 <- FindNeighbors(MALAT1, reduction = 'umap', dims = 1:2)
MALAT1 <- FindClusters(MALAT1, resolution = 0.05)
WT <- subset(MALAT1, idents = 0)
WT <- WT@assays$RNA@counts
WT <- WT[rowMeans(WT != 0) > 0.1,]

MALAT1 <- Read10X_h5('KO.h5')
MALAT1 <- scQC(MALAT1, mtThreshold = 0.05)
MALAT1 <- CreateSeuratObject(MALAT1)
MALAT1 <- NormalizeData(MALAT1)
MALAT1 <- FindVariableFeatures(MALAT1)
MALAT1 <- ScaleData(MALAT1)
MALAT1 <- RunPCA(MALAT1, verbose = FALSE)
MALAT1 <- RunUMAP(MALAT1, dims = 1:20)
MALAT1 <- FindNeighbors(MALAT1, reduction = 'umap', dims = 1:2)
MALAT1 <- FindClusters(MALAT1, resolution = 0.05)
KO <- subset(MALAT1, idents = 0)
KO <- KO@assays$RNA@counts
KO <- KO[rowMeans(KO != 0) > 0.1,]

writeMM(WT, 'WT.mtx')
writeLines(rownames(WT), 'genesWT.txt')
writeLines(colnames(WT), 'barcodesWT.txt')

writeMM(KO, 'KO.mtx')
writeLines(rownames(KO), 'genesKO.txt')
writeLines(colnames(KO), 'barcodesKO.txt')

MALAT1 <- scTenifoldKnk(WT, gKO = 'Malat1')
save(MALAT1, file = 'betaMALATko.RData')