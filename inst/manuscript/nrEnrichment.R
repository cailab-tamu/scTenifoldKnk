library(enrichR)
source('https://raw.githubusercontent.com/dosorio/utilities/master/enrichments/collapseEnrichrPathways.R')

TREM2 <- read.csv('TREM2/Results/drTREM2.csv', row.names = 1)
gList <- TREM2$gene[TREM2$p.adj < 0.05]
E <- collapseEnrichrPathways(geneList = gList,
                             dbList =  c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"),
                             jaccardThreshold = 0.99)
write.csv(E, file = 'TREM2/Results/nr_enrichmentDR_Trem2.csv', row.names = FALSE)

CFTR <- read.csv('CFTR/Results/scTenifoldKnk_SRS4245406.csv')
gList <- CFTR$gene[CFTR$p.adj < 0.05]
E <- collapseEnrichrPathways(geneList = gList,
                             dbList =  c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"),
                             jaccardThreshold = 0.99)
write.csv(E, file = 'CFTR/Results/nr_enrichmentDR_CFTR.csv', row.names = FALSE)


NKX21 <- read.csv('NKX2-1/Results/dr_GSM3716703.csv')
gList <- NKX21$gene[NKX21$p.adj < 0.05]
E <- collapseEnrichrPathways(geneList = gList,
                             dbList =  c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"),
                             jaccardThreshold = 0.99)
write.csv(E, file = 'NKX2-1/Results/nr_enrichmentDR_Nkx2-1.csv', row.names = FALSE)

HNF4AG <- read.csv('HNF4A-HNF4G/Results/drHNF4AG.csv', row.names = 1)
gList <- HNF4AG$gene[HNF4AG$p.adj < 0.05]
E <- collapseEnrichrPathways(geneList = gList,
                             dbList =  c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"),
                             jaccardThreshold = 0.99)
write.csv(E, file = 'HNF4A-HNF4G/Results/nr_enrichmentDR_Hnfg4ag.csv', row.names = FALSE)


DMD <- read.csv('DMD/Results/scTenifoldKnk_SRS4245406.csv')
gList <- DMD$gene[DMD$p.adj < 0.05]
E <- collapseEnrichrPathways(geneList = gList,
                             dbList =  c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"),
                             jaccardThreshold = 0.99)
write.csv(E, file = 'DMD/Results/nr_enrichmentDR_DMD.csv', row.names = FALSE)

MECP2r1 <- read.csv('MECP2/Results/drSRS3059998.csv', row.names = 1)
gList <- MECP2r1$gene[MECP2r1$p.adj < 0.05]
E <- collapseEnrichrPathways(geneList = gList,
                             dbList =  c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"),
                             jaccardThreshold = 0.99)
write.csv(E, file = 'MECP2/Results/nr_enrichmentDR_MECP2r1.csv', row.names = FALSE)

MECP2r2<- read.csv('MECP2/Results/drSRS3059999.csv', row.names = 1)
gList <- MECP2r2$gene[MECP2r2$p.adj < 0.05]
E <- collapseEnrichrPathways(geneList = gList,
                             dbList =  c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"),
                             jaccardThreshold = 0.99)
write.csv(E, file = 'MECP2/Results/nr_enrichmentDR_MECP2r2.csv', row.names = FALSE)
