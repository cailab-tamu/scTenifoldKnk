source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')

load('sccReal.RData')
realData <- R

load('sccSim.RData')
simData <- R

rm(R)

library(patchwork)
realList <- realData$diffRegulation$gene[realData$diffRegulation$p.adj < 0.05]
simList <- simData$diffRegulation$gene[simData$diffRegulation$p.adj < 0.05]
sharedList <- intersect(realList, simList)
length(sharedList)/length(unique(c(realList, simList)))

png(width = 3000, height = 1500, res = 300, filename = 'dr_Trem2.png')
plotDR(realData, boldGenes = sharedList, title = 'sccTenifoldNet') | plotDR(simData, boldGenes = sharedList, title = 'sccTenifoldKnk')
dev.off()

png(width = 3000, height = 3000, res = 300, filename = 'ego_Trem2.png', pointsize = 20)
plotKO(simData, 'Trem2', q = 0.9)
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
allDB <- c(BP,CC,MF)

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
