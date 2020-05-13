load('Preenterocytes.RData')
DR <- O$diffRegulation
DR <- DR[!grepl('^Rpl|^Rps|^mt\\-',DR$gene, ignore.case = TRUE),]
DR$FC <- DR$distance^2/mean(DR$distance[-1]^2)
DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
DR$p.adj <- p.adjust(DR$p.value,method = 'fdr')
drGenes <- DR$gene[DR$p.adj < 0.05]

library(Matrix)
KO <- readMM('Preenterocytes+KO.mtx')
rownames(KO) <- readLines('geneList_Preenterocytes_KO.txt')
colnames(KO) <- readLines('barcodes_Preenterocytes_KO.txt')

WT <- readMM('Preenterocytes+WT.mtx')
rownames(WT) <- readLines('geneList_Preenterocytes_WT.txt')
colnames(WT) <- readLines('barcodes_Preenterocytes_WT.txt')

library(Seurat)
KO <- CreateSeuratObject(KO, project = 'KO')
WT <- CreateSeuratObject(WT, project = 'WT')
ALL <- merge(KO,WT)

ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)
ALL <- FindVariableFeatures(ALL)
ALL <- RunPCA(ALL)
ALL <- RunUMAP(ALL,1:20)
ALL <- RunTSNE(ALL)
library(ggplot2)
png('tSNE_preE.png', width = 1200, res = 300, height = 1000)
UMAPPlot(ALL) + xlab('UMAP 1') + ylab('UMAP 2') + theme_minimal()
dev.off()

DE <- FindMarkers(ALL, ident.1 = 'WT', ident.2 = 'KO', test.use = 'MAST', min.pct = 0.05, logfc.threshold = 0)
deGenes <- DE[!grepl('^Rpl|^Rps|^mt\\-',rownames(DE), ignore.case = TRUE),]
deGenes$ID <- rownames(deGenes)
deGenes$ID[!(abs(deGenes$avg_logFC) > 0.25 & deGenes$p_val_adj < 0.05)] <- NA

library(ggplot2)
library(ggrepel)
pointColor <- densCols(cbind(deGenes$avg_logFC,-log10(deGenes$p_val)))
pointColor[!is.na(deGenes$ID)] <- 'red'
A <- ggplot(deGenes, mapping = aes(avg_logFC, -log10(p_val), label = ID)) + geom_point(col=pointColor) + theme_minimal() + theme(plot.title = element_text(size=22)) +
  geom_text_repel(aes(fontface = ifelse(deGenes$ID %in% drGenes, 2, 1))) + 
  xlab(expression(log('Fold Change'))) + ylab(expression(-log[10]('P-value'))) + xlim(c(-1,1)) + labs(title = 'Differential Expression', subtitle = 'MAST WT - KO')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
B <- plotDR(list(diffRegulation=DR), boldGenes = intersect(deGenes$ID,drGenes)) + theme_minimal() + labs(title = 'scTenifoldKnk') + theme(plot.title = element_text(size=22))
library(patchwork)
png('Ahr_preE.png', width = 3600,height = 1800, res = 300)
A + B 
dev.off()

deGenes <- deGenes$ID
intersect(deGenes,drGenes)

library(enrichR)
E <- listEnrichrDbs()
E <- enrichr(genes = deGenes, databases = c('Reactome_2016','KEGG_2019_Mouse','BioPlanet_2019','GO_Biological_Process_2018'))
E <- do.call(rbind.data.frame,E)
E <- E[E$Adjusted.P.value < 0.05,]
write.csv(E,'realAhrEnrichr_preE.csv')

E <- listEnrichrDbs()
E <- enrichr(genes = drGenes, databases = c('Reactome_2016','KEGG_2019_Mouse','BioPlanet_2019','GO_Biological_Process_2018'))
E <- do.call(rbind.data.frame,E)
E <- E[E$Adjusted.P.value < 0.05,]
write.csv(E,'simAhrEnrichr_preE.csv')

library(fgsea)
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
CHEA <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=ChEA_2016')
BIOP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
MGI <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=MGI_Mammalian_Phenotype_Level_4_2019')
PANGLAODB <- read.csv('PanglaoDB_markers_27_Mar_2020.tsv', sep = '\t', stringsAsFactors = FALSE)
cellNames <- unique(PANGLAODB$cell.type)
PANGLAODB <- sapply(cellNames, function(X){
  PANGLAODB$official.gene.symbol[PANGLAODB$cell.type %in% X]
})

Z1 <- DR$Z
names(Z1) <- toupper(DR$gene)
Z2 <- abs(DE$avg_logFC)
names(Z2) <- toupper(rownames(DE))

library(ggplot2)
set.seed(1)
T1 <- fgseaMultilevel(PANGLAODB,Z1)
set.seed(1)
T2 <- fgseaMultilevel(PANGLAODB,Z2)

png('preE_simpdb.png', height = 1000, width = 1000, res = 300)
plotEnrichment(PANGLAODB$Enterocytes,Z1) + labs(title = 'Enterocytes - scTenifoldKnk', subtitle = paste0('FDR = ',formatC(T1$padj[T1$pathway == 'Enterocytes'],digits = 2))) + xlab('Gene Rank') + ylab('Enrichment Score')
dev.off()
png('preE_realpdb.png', height = 1000, width = 1000, res = 300)
plotEnrichment(PANGLAODB$Enterocytes,Z2) + labs(title = 'Enterocytes - Abs(FC) real', subtitle = paste0('FDR = ',formatC(T2$padj[T2$pathway == 'Enterocytes'],digits = 2))) + xlab('Gene Rank') + ylab('Enrichment Score')
dev.off()

set.seed(1)
E1 <- fgsea(REACTOME, Z1,1e6)
set.seed(1)
E2 <- fgsea(REACTOME, Z2,1e6)

E1 <- E1[E1$NES > 0 & E1$padj < 0.05,]
E1$leadingEdge <- unlist(lapply(E1$leadingEdge, function(X){paste0(sort(X), collapse = ';')}))
E2 <- E2[E2$padj < 0.05,]
E2$leadingEdge <- unlist(lapply(E2$leadingEdge, function(X){paste0(sort(X), collapse = ';')}))
write.csv(E1, 'E_preE_AhrKO_sim.csv')
write.csv(E2, 'E_preE_AhrKO_real.csv')

I <- intersect(E1$pathway, E2$pathway)
writeLines(I,'I_Ahr.txt')

library(UpSetR)
png('upset_preE.png', width = 1000, height = 1500, res = 300)
upset(fromList(list(Real = E2$pathway, Simulation = E1$pathway)))
dev.off()
# gList1 <- SRS3059998$diffRegulation$gene
# gList2 <- SRS3059999$diffRegulation$gene
# 
# 
# overlapR <- pbsapply(seq_len(10), function(Z){
#   nRanks <- min(c(length(gList1), length(gList2)))
#   rList1 <- sample(gList1)
#   rList2 <- sample(gList2)
#   rOverlap <- sapply(seq_len(nRanks), function(X){
#     length(intersect(rList1[seq_len(X)], rList2[seq_len(X)]))
#   })
#   return(rOverlap)
# })
# 
# nRanks <- min(c(length(gList1), length(gList2)))
# rOverlap <- sapply(seq_len(nRanks), function(X){
#   length(intersect(gList1[seq_len(X)], gList2[seq_len(X)]))
# })