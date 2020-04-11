setwd("~/scTenifoldKNK/inst/benchmark/SERGIO")
A <- read.csv('A.csv', header = FALSE)
round(mean(A == 0),2)
corMatrix <- cor(t(A), method = 'sp')
diag(corMatrix) <- 0
corMatrix <- corMatrix/max(abs(corMatrix))
diag(corMatrix) <- 1
Matrix::image(Matrix::Matrix(corMatrix))