library(scTenifoldNet)
library(scTenifoldKnk)
load('SRS3161261_Pulmonary alveolar type II cellsAkap7.RData')
library(fgsea)
REACTOME <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')

gList <- c('Cftr','Akap7', 'Zranb1', 'Krcc1', 'Mta1', 'Rps8')

pProfiles <- lapply(gList, function(geneName){
  X <- R$WT
  Y <- X
  Y[geneName,] <- 0
  MA <- manifoldAlignment(X,Y,2)
  DR <- scTenifoldKnk:::dRegulation(MA,1)
  Z <- DR$Z
  names(Z) <- toupper(DR$gene)
  
  return(Z)
})

pProfiles <- lapply(pProfiles, function(Z){
  Z[toupper(rownames(X))]
})

tProfile <- (as.data.frame(pProfiles))
colnames(tProfile) <- gList
write.csv(tProfile, 'negControls.csv')

sGenes <- lapply(pProfiles, function(Z){
  names(Z[order(Z, decreasing = TRUE)][2:50])
})

tProfile <- t(t(tProfile)/apply(abs(tProfile),2,max))

library(ComplexHeatmap)
cftrGenes <- rownames(tProfile[order(tProfile[,1], decreasing = TRUE),][2:25,])
allGenes <- unique(unlist(sGenes))

HM <- Heatmap(tProfile[unique(unlist(sGenes)),], show_row_names = FALSE, show_row_dend = FALSE, show_heatmap_legend = FALSE, show_column_dend = FALSE) + 
  rowAnnotation(link = anno_mark(at = which(allGenes %in% cftrGenes), 
                                 labels = cftrGenes, 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
png('negHM.png', width = 800, height = 1000, res = 300)
print(HM)
dev.off()

PCA <- as.data.frame(prcomp(t(tProfile))$x)

library(ggplot2)
library(ggrepel)

png('negControls.png', width = 1100, height = 1000, res = 300)
ggplot(PCA, aes(PC1, PC2, label = rownames(PCA))) + geom_point() + geom_text_repel() + theme_bw()
dev.off()
