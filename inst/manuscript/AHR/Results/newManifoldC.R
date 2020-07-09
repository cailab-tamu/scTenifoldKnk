library(scTenifoldKnk)
library(Matrix)

load('Preenterocytes.RData')
X <- t(O$WT)
Y <- X
Y['Ahr',] <- 0
MA <- scTenifoldKnk:::manifoldAlignment(X,Y)
MA <- MA[!grepl('_RPL|_RPS|_RP[[:digit:]]+',rownames(MA), ignore.case = TRUE),]
DR <- scTenifoldKnk:::dRegulation(MA, 'Ahr')
O$manifoldAlignment <- MA
O$diffRegulation <- DR
O$WT <- X
O$KO <- Y
save(O, file = 'YC_Preenterocytes.RData')
