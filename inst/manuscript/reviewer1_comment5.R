library(enrichR)

dumpEnrichment <- function(genesDR){
  E <- do.call(rbind.data.frame, enrichr(genesDR, databases = c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")))
  E <- E[E$Adjusted.P.value < 0.05,]
  E <- E[order(E$Odds.Ratio, decreasing = TRUE),]
  rownames(E) <- NULL
  return(E)
}

# Trem2
DR <- read.csv('TREM2/Results/drTREM2.csv', row.names = 1)
genesDR <- DR$gene[DR$p.adj < 0.05]
write.csv(dumpEnrichment(genesDR), 'TREM2/Results/enrichmentDR_Trem2.csv')

# Nkx2-1
DR <- read.csv('NKX2-1/Results/dr_GSM3716703.csv', row.names = 1)
genesDR <- DR$gene[DR$p.adj < 0.05]
write.csv(dumpEnrichment(genesDR), 'NKX2-1/Results/enrichmentDR_Nkx2-1.csv')

# Hnfg4ag
DR <- read.csv('HNF4A-HNF4G/Results/drHNF4AG.csv', row.names = 1)
genesDR <- DR$gene[DR$p.adj < 0.05]
write.csv(dumpEnrichment(genesDR), 'HNF4A-HNF4G/Results/enrichmentDR_Hnfg4ag.csv')

# Cftr
DR <- read.csv('CFTR/Results/scTenifoldKnk_SRS4245406.csv')
genesDR <- DR$gene[DR$p.adj < 0.05]
write.csv(dumpEnrichment(genesDR), 'CFTR/Results/enrichmentDR_CFTR.csv')

# Dmd
DR <- read.csv('DMD/Results/scTenifoldKnk_SRS4245406.csv')
genesDR <- DR$gene[DR$p.adj < 0.05]
write.csv(dumpEnrichment(genesDR), 'DMD/Results/enrichmentDR_DMD.csv')

# Mecp2
DR <- read.csv('MECP2/Results/drSRS3059998.csv')
genesDR <- DR$gene[DR$p.adj < 0.05]
write.csv(dumpEnrichment(genesDR), 'MECP2/Results/enrichmentDR_MECP2r1.csv')

DR <- read.csv('MECP2/Results/drSRS3059999.csv')
genesDR <- DR$gene[DR$p.adj < 0.05]
write.csv(dumpEnrichment(genesDR), 'MECP2/Results/enrichmentDR_MECP2r2.csv')
