source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')

load('sccReal.RData')
realData <- R

load('sccSim.RData')
simData <- R

rm(R)

realList <- realData$diffRegulation$gene[realData$diffRegulation$p.adj < 0.05]
simList <- simData$diffRegulation$gene[simData$diffRegulation$p.adj < 0.05]
sharedList <- intersect(realList, simList)

realOrder <- realData$diffRegulation$gene
simOrder <- simData$diffRegulation$gene

bootOrder <- pbapply::pbsapply(1:10, function(X){
  R <- sample(realOrder)
  S <- sample(simOrder)
  sapply(seq_len(min(c(length(realOrder), length(simOrder)))), function(i){
  length(intersect(R[seq_len(i)], S[seq_len(i)]))
  })
})

ranks <- seq_len(min(c(length(realOrder), length(simOrder))))
overlaps <- sapply(ranks, function(i){
  length(intersect(realOrder[seq_len(i)], simOrder[seq_len(i)]))
})

oList <- data.frame(Rank = ranks, Overlap = overlaps)
rList <- data.frame(Rank = ranks, Overlap = rowMeans(bootOrder), Min = apply(bootOrder,1,min), Max = apply(bootOrder,1,max))

png('oL_Trem2.png', width = 1000, height = 1200, res = 300, pointsize = 20)
library(ggplot2)
ggplot(oList, aes(Rank, Overlap)) +
  geom_line(aes(color='blue')) +
  theme_bw() + ylab('Size of Overlap') +
  geom_line(data = rList, aes(color = 'orange'), lty = 2) +
  geom_ribbon(data = rList, mapping = aes(ymin = Min, ymax = Max), alpha = 0.3, fill = 'orange') +
  scale_color_manual(values=c("#999999", "#E69F00"),
                    name="Order\nComparison",
                    breaks=c("blue", "orange"),
                    labels=c("Real", "Random")) + theme(legend.position="bottom")
dev.off()

png(width = 3000, height = 1500, res = 300, filename = 'dr_Trem2.png')
library(patchwork)
plotDR(realData, boldGenes = sharedList, title = 'sccTenifoldNet') | plotDR(simData, boldGenes = sharedList, title = 'sccTenifoldKnk')
dev.off()

png(width = 3000, height = 3000, res = 300, filename = 'ego_Trem2.png', pointsize = 20)
plotKO(simData, 'Trem2', q = 0.99)
dev.off()

library(UpSetR)
png(width = 750, height = 750, res = 300, filename = 'gc_Trem2.png')
upset(fromList(list(Real=realList, Simulated = simList)))
dev.off()


library(fgsea)
BP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018')
CC <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2018')
MF <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Molecular_Function_2018')
KEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
BIOP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
REAC <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
allDB <- c(BP,CC,MF, KEGG, BIOP, REAC)

zValues <- simData$diffRegulation$Z
names(zValues) <- toupper(simData$diffRegulation$gene)
gseaResult <- fgseaMultilevel(allDB, zValues)
gseaResult <- gseaResult[gseaResult$padj < 0.05 & gseaResult$NES > 0,]
gseaResult <- gseaResult[order(gseaResult$NES, decreasing = TRUE),]

png(width = 1000, height = 1000, res = 300, filename = 'GSEA1.png')
P <- which(gseaResult$pathway %in% 'regulation of cholesterol metabolic process (GO:0090181)')
fgsea::plotEnrichment(pathway = allDB[['regulation of cholesterol metabolic process (GO:0090181)']], zValues) +
  ggtitle('Regulation of Cholesterol \nMetabolic Process (GO:0090181)', subtitle = paste0('P = ', formatC(gseaResult$pval[P], digits = 3, format = 'e'))) +
  xlab('Gene Rank') + ylab('Enrichment Score')
dev.off()

png(width = 1000, height = 1000, res = 300, filename = 'GSEA2.png')
P <- which(gseaResult$pathway %in% 'autolysosome (GO:0044754)')
fgsea::plotEnrichment(pathway = allDB[['autolysosome (GO:0044754)']], zValues) +
  ggtitle('Autolysosome\n(GO:0044754)', subtitle = paste0('P = ', formatC(gseaResult$pval[1], digits = 3, format = 'e'))) +
  xlab('Gene Rank') + ylab('Enrichment Score')
dev.off()

png(width = 1000, height = 1000, res = 300, filename = 'GSEA3.png')
P <- which(gseaResult$pathway %in% 'chylomicron assembly (GO:0034378)')
fgsea::plotEnrichment(pathway = allDB[['chylomicron assembly (GO:0034378)']], zValues) +
  ggtitle('Chylomicron Assembly\n(GO:0034378)', subtitle = paste0('P = ', formatC(gseaResult$pval[1], digits = 3, format = 'e'))) +
  xlab('Gene Rank') + ylab('Enrichment Score')
dev.off()


realZ <- realData$diffRegulation$distance
names(realZ) <- realData$diffRegulation$gene

simZ <- simData$diffRegulation$distance #simData$diffRegulation$Z
names(simZ) <- simData$diffRegulation$gene

sharedGenes <- intersect(realData$diffRegulation$gene, simData$diffRegulation$gene)
round(cor(realZ[sharedGenes], simZ[sharedGenes], method = 'sp'),2)

Z <- data.frame(realZ = realZ[sharedGenes], simZ = simZ[sharedGenes], geneName = sharedGenes)
Z$geneName[!Z$geneName %in% sharedList] <- NA
pColor <-densCols(Z[,1:2])
pColor[!is.na(Z$geneName)] <- 'red'
library(ggplot2)
library(ggrepel)
library(statsExpressions)
png('cor_Trem2.png', width = 1500, height = 1500, res = 300, pointsize = 20)
ggplot(Z, aes(realZ, simZ, label = geneName)) +
  geom_point(col= pColor) +
  theme_bw() +
  geom_text_repel() +
  ylab('Z-score Simulation') +
  xlab('Z-score Real') +
  labs(subtitle = expression(paste(NULL, "log"["e"](italic("S")), " = ",
                                   "20.76", ", ", italic("p"), " = ",
                                   "< 0.001", ", ", widehat(italic(rho))["Spearman"],
                                   " = ", "0.79")))#expr_corr_test(data = Z, x = realZ, y = simZ, type = 'spearman'))
dev.off()
