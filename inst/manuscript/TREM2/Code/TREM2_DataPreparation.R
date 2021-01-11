library(Matrix)

# Reading input files
cellInfo <- read.csv('../Data/GSE130626_cell_info.csv.gz', stringsAsFactors = FALSE)
geneInfo <- read.csv('../Data/GSE130626_gene_info.csv.gz', stringsAsFactors = FALSE)
geneInfo <- geneInfo[complete.cases(geneInfo),]
UMI <- read.csv('../Data/GSE130626_umi_counts.csv.gz', row.names = 1)

# Replacing Ensembl to Symbols
UMI <- Matrix(as.matrix(UMI))
UMI <- UMI[geneInfo$gene_id,]
rownames(UMI) <- geneInfo$symbol

# Selecting treatment
cellInfo <- cellInfo[cellInfo$treatment == 'cuprizone',]

# Selecting cells
WT <- cellInfo$cell_id[cellInfo$trem2_genotype == 'WT']
KO <- cellInfo$cell_id[cellInfo$trem2_genotype == 'KO']

# Filtering cells
WT <- UMI[,WT]
KO <- UMI[,KO]

# Writing objects
writeMM(WT, file = '../Data/matrix_WT.mtx')
writeMM(KO, file = '../Data/matrix_KO.mtx')

writeLines(rownames(WT), '../Data/geneNames_WT.txt')
writeLines(rownames(KO), '../Data/geneNames_KO.txt')

writeLines(colnames(WT), '../Data/barcodes_WT.txt')
writeLines(colnames(KO), '../Data/barcodes_KO.txt')
