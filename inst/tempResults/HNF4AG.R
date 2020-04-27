fixPValues <- function(X){
  X$manifoldAlignment <- X$manifoldAlignment[,1:2]
  X$diffRegulation <- scTenifoldNet::dRegulation(X$manifoldAlignment)
  X$diffRegulation <- X$diffRegulation[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^mt-', X$diffRegulation$gene, ignore.case = TRUE),]
  D <- X$diffRegulation$distance[-c(1:2)]
  X$diffRegulation$FC <- X$diffRegulation$distance^2/mean(D^2)
  X$diffRegulation$p.value <- pchisq(X$diffRegulation$FC, df = 1, lower.tail = FALSE)
  X$diffRegulation$p.adj <- p.adjust(X$diffRegulation$p.value, method = 'fdr')
  return(X)
}
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')

load('all0.9/HNF4AG_0.9_GSM4116571.RData')
GSM3477499 <- fixPValues(GSM3477499)
save(GSM3477499, file = '../manuscript/HNF4A-HNF4G/Results/GSM3477499.RData')

plotDR(GSM3477499)
