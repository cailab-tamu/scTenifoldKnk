library(Matrix)
library(Seurat)
library(harmony)

WT <- readMM('../Data/Preenterocytes_WT.mtx')
rownames(WT) <- readLines('../Data/Preenterocytes_WT_genes.txt')
colnames(WT) <- readLines('../Data/Preenterocytes_WT_barcodes.txt')
WT <- WT[!grepl('^MT-|^RPL|^RPS',rownames(WT), ignore.case = TRUE),]

KO <- readMM('../Data/Preenterocytes_KO.mtx')
rownames(KO) <- readLines('../Data/Preenterocytes_KO_genes.txt')
colnames(KO) <- readLines('../Data/Preenterocytes_KO_barcodes.txt')
KO <- KO[!grepl('^MT-|^RPL|^RPS',rownames(KO), ignore.case = TRUE),]

WT <- CreateSeuratObject(WT, project = 'WT')
KO <- CreateSeuratObject(KO, project = 'KO')

ALL <- merge(WT,KO)
ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)
ALL <- FindVariableFeatures(ALL)
ALL <- RunPCA(ALL)
ALL <- RunHarmony(ALL, group.by.vars = 'orig.ident')
ALL <- RunTSNE(ALL, reduction = 'harmony')
ALL <- RunUMAP(ALL, reduction = 'harmony', dims = 1:50)
UMAPPlot(ALL)

DE <- FindMarkers(ALL, ident.1 = 'KO', ident.2 = 'WT', test.use = 'MAST', logfc.threshold = 0)
DE <- DE[order(abs(DE$avg_logFC), decreasing = TRUE),]
deList <- rownames(DE)
deZ <- abs(DE$avg_logFC)
names(deZ) <- toupper(rownames(DE))

load('../Results/Preenterocytes.RData')
# MA <- O$manifoldAlignment
# MA <- MA[!grepl('_MT-|_RPL|_RPS',rownames(MA), ignore.case = TRUE),]
# DR <- scTenifoldNet::dRegulation(MA)
DR$FC <- (DR$distance^2)/mean(DR$distance[-1]^2)
DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
O$diffRegulation <- DR
drZ <- DR$Z
names(drZ) <- toupper(DR$gene)

drList <- DR$gene
deList <- deList[deList %in% drList]
drList <- drList[drList %in% deList]

library(OrderedList)
png('DE_DR.png', width = 2000, height = 1000, res = 300)
par(mar=c(3,3,1,1), mgp=c(1.5,0.5,0))
plot(compareLists(deList,drList, alphas = 0.005))
dev.off()

library(fgsea)
KEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
drE <- fgsea(REACTOME, drZ, 1e6)
deE <- fgsea(REACTOME, deZ, 1e6)
library(UpSetR)
png('fgsea.png', width = 600, height = 600, res = 300)
upset(fromList(list(Simulation=drE$pathway[drE$padj < 0.05 & drE$NES > 0],Real=deE$pathway[deE$padj < 0.05 & deE$NES > 0])))
dev.off()

library(ggplot2)
library(ggrepel)
pData <- data.frame(FC=DE$avg_logFC, P=-log10(DE$p_val), G=rownames(DE))
pData$G[!pData$G %in% DR$gene[DR$p.adj < 0.05]] <- NA
pData$L <- seq_len(nrow(pData))
pData$L[!is.na(pData$G)] <- 1e6
pData <- pData[order(pData$L, decreasing = FALSE),]
pData$C <- densCols(pData[,1:2])
pData$C[!is.na(pData$G)] <- 'red'
png('DE_drLabels.png', width = 2000, height = 2000, res = 300)
ggplot(pData, aes(FC,P,label=G)) + geom_point(color=pData$C) + geom_text_repel() + theme_bw() + ylim(c(0,80)) + xlim(c(-.6,.6))
dev.off()       

source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')

png('KO.png', width = 3500, height = 3500, res = 300)
plotKO(O, 'Ahr', nCategories = 10)
dev.off()
