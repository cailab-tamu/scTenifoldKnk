#' @importFrom Matrix Matrix colSums
#' @importFrom stats lm predict
scQC <- function(X, mtThreshold = 0.1){
  if(class(X) == 'Seurat'){
    countMatrix <- X@assays$RNA@counts
  } else {
    countMatrix <- X
  }
  librarySize <- colSums(countMatrix)
  mtCounts <- colSums(countMatrix[grep('^MT-',toupper(rownames(countMatrix))),])
  nGenes <- colSums(countMatrix != 0)
  mtProportion <- mtCounts/librarySize

  genesLM <- lm(nGenes~librarySize)
  mtLM <- lm(mtCounts~librarySize)

  genesLM <- as.data.frame(predict(genesLM, data.frame(librarySize), interval = 'prediction'))
  mtLM <- as.data.frame(predict(mtLM, data.frame(librarySize), interval = 'prediction'))

  selectedCells <- mtCounts > mtLM$lwr & mtCounts < mtLM$upr & nGenes > genesLM$lwr & nGenes < genesLM$upr & mtProportion <= mtThreshold & librarySize < 2 * mean(librarySize)
  selectedCells <- colnames(countMatrix)[selectedCells]
  if(class(X) == 'Seurat'){
    X <- subset(X, cells = selectedCells)
  } else {
    X <- countMatrix[,selectedCells]
  }
  return(X)
}
