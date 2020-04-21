fixPValues <- function(X){
  X$manifoldAlignment <- X$manifoldAlignment[,1:2]
  X$diffRegulation <- scTenifoldNet::dRegulation(X$manifoldAlignment)
  X$diffRegulation <- X$diffRegulation[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^mt-', X$diffRegulation$gene, ignore.case = TRUE),]
  D <- X$diffRegulation$distance[-1]
  X$diffRegulation$FC <- X$diffRegulation$distance^2/mean(D^2)
  X$diffRegulation$p.value <- pchisq(X$diffRegulation$FC, df = 1, lower.tail = FALSE)
  X$diffRegulation$p.adj <- p.adjust(X$diffRegulation$p.value)
  return(X)
}
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')

load('all0.9/MECP2_0.9_SRS3059998.RData')
SRS3059998 <- fixPValues(WT)
save(SRS3059998, file = '../manuscript/MECP2/SRS3059998.RData')

load('all0.9/MECP2_0.9_SRS3059999.RData')
SRS3059999 <- fixPValues(WT)
save(SRS3059999, file = '../manuscript/MECP2/SRS3059999.RData')
