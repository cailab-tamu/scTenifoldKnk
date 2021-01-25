library(Matrix)
countMatrix <- read.csv('simulationOutput.csv', header = FALSE)
countMatrix <- as.matrix(countMatrix)
countMatrix <- Matrix(countMatrix)
mean(countMatrix != 0)

rownames(countMatrix) <- paste0('g', seq_len(nrow(countMatrix)))
colnames(countMatrix) <- paste0('c', seq_len(ncol(countMatrix)))

set.seed(1)
O <- scTenifoldKnk(countMatrix, gKO = 20, nc_nComp = 3, nc_q = 0.85, qc_mtThreshold = 1, qc_minLSize = 0, nc_lambda = 0.7)
DR <- O$diffRegulation
plot(DR[paste0(1:100),]$Z)
