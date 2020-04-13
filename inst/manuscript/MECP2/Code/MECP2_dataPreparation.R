library(Matrix)
library(XML)
library(xml2)

cellID <- read.csv('../Data/SRA667466_SRS3059998.clusters.txt', header = FALSE, sep = ' ')
download_html('https://panglaodb.se/list_clusters_and_cell_types.html?sra=SRA667466&srs=SRS3059998', file = '../Data/SRA667466_SRS3059998.clusters')
clusterID <- readHTMLTable('../Data/SRA667466_SRS3059998.clusters', as.data.frame = TRUE)
selClusters <- clusterID$tbl$`Cluster ID`[clusterID$tbl$`Inferred cell type` %in% 'Neurons']
selClusters <- as.vector(selClusters)
selCells <- as.vector(cellID$V1[cellID$V2 %in% selClusters])

load('../Data/SRA667466_SRS3059998.sparse.RData')
sm <- sm[,selCells]
rownames(sm) <- unlist(lapply(strsplit(rownames(sm), '_EN'), function(X){X[1]}))

writeMM(sm, file = '../Data/matrix_SRS3059998.mtx')
writeLines(colnames(sm), '../Data/barcodes_SRS3059998.txt')
writeLines(rownames(sm), '../Data/geneNames_SRS3059998.txt')



cellID <- read.csv('../Data/SRA667466_SRS3059999.clusters.txt', header = FALSE, sep = ' ')
download_html('https://panglaodb.se/list_clusters_and_cell_types.html?sra=SRA667466&srs=SRS3059999', file = '../Data/SRA667466_SRS3059999.clusters')
clusterID <- readHTMLTable('../Data/SRA667466_SRS3059999.clusters', as.data.frame = TRUE)
selClusters <- clusterID$tbl$`Cluster ID`[clusterID$tbl$`Inferred cell type` %in% 'Neurons']
selClusters <- as.vector(selClusters)
selCells <- as.vector(cellID$V1[cellID$V2 %in% selClusters])

load('../Data/SRA667466_SRS3059999.sparse.RData')
sm <- sm[,selCells]
rownames(sm) <- unlist(lapply(strsplit(rownames(sm), '_EN'), function(X){X[1]}))

writeMM(sm, file = '../Data/matrix_SRS3059999.mtx')
writeLines(colnames(sm), '../Data/barcodes_SRS3059999.txt')
writeLines(rownames(sm), '../Data/geneNames_SRS3059999.txt')
