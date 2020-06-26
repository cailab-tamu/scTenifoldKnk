#' @importFrom stats dist p.adjust pchisq
#' @importFrom MASS boxcox
dRegulation <- function(manifoldOutput, gKO){
  
  geneList <- rownames(manifoldOutput)
  geneList <- geneList[grepl('^X_', geneList)]
  geneList <- gsub('^X_','', geneList)
  nGenes <- length(geneList)
  
  eGenes <- nrow(manifoldOutput)/2
  
  eGeneList <- rownames(manifoldOutput)
  eGeneList <- eGeneList[grepl('^Y_', eGeneList)]
  eGeneList <- gsub('^Y_','', eGeneList)
  
  if(nGenes != eGenes){
    stop('Number of identified and expected genes are not the same')
  }
  if(!all(eGeneList == geneList)){
    stop('Genes are not ordered as expected. X_ genes should be followed by Y_ genes in the same order')
  }
  
  dMetric <- sapply(seq_len(nGenes), function(G){
    X <- manifoldOutput[G,]
    Y <- manifoldOutput[(G+nGenes),]
    I <- rbind(X,Y)
    O <- stats::dist(I)
    O <- as.numeric(O)
    return(O)
  })
  
  dMetric <- dMetric/min(dMetric)
  ### BOX-COX
  lambdaValues <- seq(-2,2,length.out = 1000)
  lambdaValues <- lambdaValues[lambdaValues != 0]
  BC <- MASS::boxcox(dMetric~1, plot=FALSE, lambda = lambdaValues)
  BC <- BC$x[which.max(BC$y)]
  if(BC < 0){
    nD <- 1/(dMetric ^ BC)
  } else {
    nD <- dMetric ^ BC
  }
  Z <- scale(nD)
  FC <- (dMetric^2)/mean((dMetric[-seq_along(gKO)]^2))
  pValues <- pchisq(q = FC,df = 1,lower.tail = FALSE)
  pAdjusted <- p.adjust(pValues, method = 'fdr')
  dOut <- data.frame(
    gene = geneList, 
    distance = dMetric,
    Z = Z,
    FC = FC,
    p.value = pValues,
    p.adj = pAdjusted
  )
  dOut <- dOut[order(dOut$p.value),]
  dOut <- as.data.frame.array(dOut)
  return(dOut)
}
