library(XML)
library(xml2)
library(fgsea)

download_html('http://amp.pharm.mssm.edu/Harmonizome/gene_set/Cystic+Fibrosis/HuGE+Navigator+Gene-Phenotype+Associations', file = 'CysticFibrosis.html')
CF <- readHTMLTable('CysticFibrosis.html')[[2]][-1,]
writeLines(as.vector(CF$Symbol), 'gSet_CF.txt')

WP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2019_Mouse')
writeLines(WP$`Mecp2 and Associated Rett Syndrome WP2910`, 'gSet_Rett.txt')


download_html('http://amp.pharm.mssm.edu/Harmonizome/gene_set/Duchenne+muscular+dystrophy+%28DMD%29_Muscle+-+Striated+%28Skeletal%29+-+Diaphragm+%28MMHCC%29_GSE1026/GEO+Signatures+of+Differentially+Expressed+Genes+for+Diseases', 'DMD.html')
DMD <- readHTMLTable('DMD.html')
DMD <- do.call(rbind.data.frame, DMD[-1])[-c(1,301),]
writeLines(unique(as.vector(DMD$Symbol)), 'gSet_DMD.txt')

