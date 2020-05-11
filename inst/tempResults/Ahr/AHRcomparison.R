load('Goblet cells.RData')
DR <- O$diffRegulation
DR <- DR[!grepl('^Rpl|^Rps|^mt\\-',DR$gene, ignore.case = TRUE),]
DR$FC <- DR$distance^2/mean(DR$distance[-1]^2)
DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
DR$p.adj <- p.adjust(DR$p.value,method = 'fdr')
drGenes <- DR$gene[DR$p.adj < 0.05]

deGenes <- read.csv('DE_Goblet cells.csv', stringsAsFactors = FALSE)
deGenes <- deGenes[!grepl('^Rpl|^Rps|^mt\\-',deGenes$X, ignore.case = TRUE),]
deGenes$ID <- deGenes$X
deGenes$ID[!(abs(deGenes$avg_logFC) > 0.15 & deGenes$p_val_adj < 0.05)] <- NA
library(ggplot2)
pointColor <- densCols(cbind(deGenes$avg_logFC,-log10(deGenes$p_val)))
pointColor[!is.na(deGenes$ID)] <- 'red'
A <- ggplot(deGenes, mapping = aes(avg_logFC, -log10(p_val), label = ID)) + geom_point(col=pointColor) + theme_minimal() + 
  geom_text_repel(aes(fontface = ifelse(deGenes$ID %in% drGenes, 2, 1))) + 
  xlab(expression(log('Fold Change'))) + ylab(expression(-log[10]('P-value'))) + xlim(c(-1,1))
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
B <- plotDR(list(diffRegulation=DR), boldGenes = intersect(deGenes$ID,drGenes)) + theme_minimal()
library(patchwork)
png('Ahr.png', width = 4000,height = 2000, res = 300)
A + B
dev.off()

deGenes <- deGenes$X[deGenes$p_val_adj < 0.05 & abs(deGenes$avg_logFC) > 0.15]
intersect(deGenes,drGenes)

library(enrichR)
E <- listEnrichrDbs()
E <- enrichr(genes = deGenes, databases = c('Reactome_2016','KEGG_2019_Mouse','BioPlanet_2019','GO_Biological_Process_2018'))
E <- do.call(rbind.data.frame,E)
E <- E[E$Adjusted.P.value < 0.05,]

library(fgsea)
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
CHEA <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=ChEA_2016')
BIOP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
MGI <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=MGI_Mammalian_Phenotype_Level_4_2019')

Z1 <- DR$Z
names(Z1) <- toupper(DR$gene)
Z2 <- deGenes$avg_logFC
names(Z2) <- toupper(deGenes$X)

set.seed(1)
E1 <- fgseaMultilevel(REACTOME, Z1)
set.seed(1)
E2 <- fgseaMultilevel(REACTOME, Z2)

E1 <- E1[E1$NES > 0 & E1$padj < 0.05,]
E2 <- E2[E2$padj < 0.05,]

intersect(E1$pathway,E2$pathway)

library(UpSetR)
png('upset.png', width = 1000, height = 1500, res = 300)
upset(fromList(list(Simulation=E1$pathway,Real=E2$pathway)))
dev.off()

# library(OrderedList)
# z1Names <- names(Z[order(Z, decreasing = TRUE)])
# z2Names <- names(Z2[order(Z2, decreasing = FALSE)])
# z1Names <- z1Names[z1Names %in% intersect(z1Names,z2Names)]
# z2Names <- z2Names[z2Names %in% intersect(z1Names,z2Names)]
# 
# plot(compareLists(z1Names,z2Names))
