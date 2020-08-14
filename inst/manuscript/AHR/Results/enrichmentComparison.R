library(fgsea)
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')

load('Preenterocytes.RData')
Z <- O$diffRegulation$Z
names(Z) <- toupper(O$diffRegulation$gene)

E1 <- fgsea(REACTOME, Z, 10)
E1 <- E1[E1$padj < 0.05,]

load('YC_Preenterocytes.RData')
Z <- O$diffRegulation$Z
names(Z) <- toupper(O$diffRegulation$gene)

E2 <- fgsea(REACTOME, Z, 10)
E2 <- E2[E2$padj < 0.05,]