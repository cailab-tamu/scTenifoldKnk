fileList <- list.files(pattern = 'trem2S')
fileList <- fileList[grepl('.csv', fileList)]

allZ <- lapply(fileList, function(X){
  fileContent <- read.csv(X)
  fileContent <- fileContent[!grepl('^Rpl|^Rps|^Rp[[:digit:]]|^Mt-', fileContent$gene),]
  Z <- fileContent$Z
  names(Z) <- fileContent$gene
  return(Z)
})

geneList <- unique(unlist(lapply(allZ, names)))

allZ <- sapply(allZ, function(X){X[geneList]})
allZ <- allZ[,allZ['Trem2',] > 0]
allZ <- allZ[complete.cases(allZ),]

library(corrplot)
corV <- cor(allZ, method = 'sp')
png('trem2Stability_SpearmanCor.png', width = 1500, height = 1500, res = 300)
corrplot.mixed(corV, tl.pos = 'd')
dev.off()

# mean(corV[lower.tri(corV)])
# sd(corV[lower.tri(corV)])

library(GSVA)
library(fgsea)

KEGG <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
rownames(allZ) <- toupper(rownames(allZ))

# E <- apply(allZ, 2, function(X){fgseaMultilevel(KEGG,X)})
# 
# cmP <- vector()
# for(i in 1:10){
#   cmP[i] <- E[[3]]$padj[E[[3]]$pathway %in% 'Cholesterol metabolism']
# }
# 
# allE <- lapply(E, function(X){
#   X$pathway[X$padj < 0.05]
# })
# 
# jaccardValues <- sapply(allE, function(i){
#   lapply(allE, function(j){
#     mean(table(c(i,j)) == 2)
#   })
# })
# 
# mean(unlist(jaccardValues[lower.tri(jaccardValues)]))
# sd(unlist(jaccardValues[lower.tri(jaccardValues)]))
# 
# 
# rownames(allZ) <- toupper(rownames(allZ))
E <- gsva(allZ, KEGG, method='ssgsea', min.sz=2)

library(ComplexHeatmap)
png('trem2enrichmentHeatmap.png', height = 2000)
Heatmap(E)
dev.off()

R <- sapply(1:100, function(X){
  rownames(allZ) <- sample(rownames(allZ))
  E <- gsva(allZ, KEGG, method='ssgsea')
  E <- E['Cholesterol metabolism',]
  return(E)
})

R <- t(R)
colnames(R) <- paste0('S', seq_len(ncol(R)))
R <- as.data.frame(R)

library(ggplot2)
ggplot(R, aes(S2)) + geom_density() + geom_vline(xintercept = E['Cholesterol metabolism',2], lty = 2, col = 'red') + theme_bw() + ylab('Density')


# # allG <- lapply(fileList, function(X){
# #   fileContent <- read.csv(X)
# #   fileContent <- fileContent[!grepl('^Rpl|^Rps|Mt-|^Rp[[:digit:]]',fileContent$gene),]
# #   fileContent$FC <- fileContent$distance^2/mean(fileContent$distance[-1]^2)
# #   fileContent$p.value <- pchisq(fileContent$FC, df = 1, lower.tail = FALSE)
# #   fileContent$p.adj <- p.adjust(fileContent$p.value, method = 'fdr')
# #   G <- fileContent$gene[fileContent$p.adj < 0.05]
# #   return(G)
# # })
# # 
# # for(i in seq_along(allG)){
# #   for(j in seq_along(allG)){
# #     gList <- c(allG[[i]], allG[[j]])
# #     print(mean(table(gList) == 2))
# #   }
# # }
