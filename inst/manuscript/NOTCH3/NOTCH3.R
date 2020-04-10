library(Matrix)
GSM4314258 <- readMM('GSM4314258_counts_sample1.mtx.gz')
colnames(GSM4314258) <- readLines('GSM4314258_barcodes_sample1.txt.gz')
rownames(GSM4314258) <- readLines('GSM4314258_genes_sample1.txt.gz')

GSM4314259 <- readMM('GSM4314259_counts_sample2.mtx.gz')
colnames(GSM4314259) <- readLines('GSM4314259_barcodes_sample2.txt.gz')
rownames(GSM4314259) <- readLines('GSM4314259_genes_sample2.txt.gz')

