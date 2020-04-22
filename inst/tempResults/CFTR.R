fixPValues <- function(X){
  X$manifoldAlignment <- X$manifoldAlignment[,1:2]
  X$diffRegulation <- scTenifoldNet::dRegulation(X$manifoldAlignment)
  X$diffRegulation <- X$diffRegulation[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^mt-', X$diffRegulation$gene, ignore.case = TRUE),]
  D <- X$diffRegulation$distance[-1]
  X$diffRegulation$FC <- X$diffRegulation$distance^2/mean(D^2)
  X$diffRegulation$p.value <- pchisq(X$diffRegulation$FC, df = 1, lower.tail = FALSE)
  X$diffRegulation$p.adj <- p.adjust(X$diffRegulation$p.value, method = 'fdr')
  return(X)
}
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')


load('all0.9/CFTR_0.9_SRS3161261.RData')
SRS3161261 <- fixPValues(SRS3161261)
save(SRS3161261, file = '../manuscript/CFTR/Results/SRS3161261.RData')

load('all0.9/CFTR_0.9_SRS4245406.RData')
SRS4245406 <- fixPValues(SRS4245406)
save(SRS4245406, file = '../manuscript/CFTR/Results/SRS4245406.RData')

