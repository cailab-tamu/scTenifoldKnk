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

load('all0.9/NKX21_0.9_GSM3716703.RData')
GSM3716703 <- fixPValues(WT)
save(GSM3716703, file = '../manuscript/NKX2-1/Results/GSM3716703.RData')

plotDR(GSM3716703)

