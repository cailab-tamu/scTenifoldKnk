library(Matrix)
library(XML)
library(xml2)
cellID <- read.csv('../Data/SRA688150_SRS3161261.clusters.txt', header = FALSE, sep = ' ')
download_html('https://panglaodb.se/list_clusters_and_cell_types.html?sra=SRA688150&srs=SRS3161261', file = '../Data/SRA688150_SRS3161261.clusters')
clusterID <- readHTMLTable('../Data/SRA688150_SRS3161261.clusters', as.data.frame = TRUE)
selClusters <- clusterID$tbl$`Cluster ID`[clusterID$tbl$`Inferred cell type` %in% 'Pulmonary alveolar type II cells']
selClusters <- as.vector(selClusters)
selCells <- as.vector(cellID$V1[cellID$V2 %in% selClusters])

load('../Data/SRA688150_SRS3161261.sparse.RData')
sm <- sm[,selCells]
rownames(sm) <- unlist(lapply(strsplit(rownames(sm), '_EN'), function(X){X[1]}))

writeMM(sm, file = '../Data/matrix_SRS3161261.mtx')
writeLines(colnames(sm), '../Data/barcodes_SRS3161261.txt')
writeLines(rownames(sm), '../Data/geneNames_SRS3161261.txt')


cellID <- read.csv('../Data/SRA692557_SRS4245406.clusters.txt', header = FALSE, sep = ' ')
download_html('https://panglaodb.se/list_clusters_and_cell_types.html?sra=SRA692557&srs=SRS4245406', file = '../Data/SRA692557_SRS4245406.clusters')
clusterID <- readHTMLTable('../Data/SRA692557_SRS4245406.clusters', as.data.frame = TRUE)
selClusters <- clusterID$tbl$`Cluster ID`[clusterID$tbl$`Inferred cell type` %in% 'Pulmonary alveolar type II cells']
selClusters <- as.vector(selClusters)
selCells <- as.vector(cellID$V1[cellID$V2 %in% selClusters])

load('../Data/SRA692557_SRS4245406.sparse.RData')
sm <- sm[,selCells]
rownames(sm) <- unlist(lapply(strsplit(rownames(sm), '_EN'), function(X){X[1]}))

writeMM(sm, file = '../Data/matrix_SRS4245406.mtx')
writeLines(colnames(sm), '../Data/barcodes_SRS4245406.txt')
writeLines(rownames(sm), '../Data/geneNames_SRS4245406.txt')
