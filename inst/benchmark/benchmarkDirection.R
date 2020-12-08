library(scTenifoldNet)
library(ggplot2)

SERGIO <- read.csv('SERGIO/simulationOutput.csv', header = FALSE)
rownames(SERGIO) <- paste0('G', seq_len(nrow(SERGIO)))
colnames(SERGIO) <- paste0('C', seq_len(ncol(SERGIO)))

dirDiff <- sapply(2:99, function(V){
  dirMatrix <- pcNet(as.matrix(SERGIO), nComp = V)

  symMatrix <- (dirMatrix + t(dirMatrix))/2
  mean((dirMatrix - symMatrix)^2)
})

DF <- data.frame(nPC = (seq_along(dirDiff)+1), dirDiff = dirDiff)
png('dirDifference.png',width = 1000, height = 1000, res= 300)
ggplot(DF, aes(nPC, dirDiff)) + 
  geom_point() + 
  geom_smooth() + 
  theme_bw() + 
  ylab(Average~Difference~from~Symmetric) + 
  xlab(Number~of~PCs)
dev.off()

dirDiff <- sapply(2:99, function(V){
  set.seed(1)
  dirMatrix <- makeNetworks(as.matrix(SERGIO), nComp = V, q = 0)
  set.seed(1)
  dirMatrix <- tensorDecomposition(dirMatrix, nDecimal = 3)$X
  
  symMatrix <- (dirMatrix + t(dirMatrix))/2
  mean((dirMatrix - symMatrix)^2)
})

DF <- data.frame(nPC = (seq_along(dirDiff)+1), dirDiff = dirDiff)
png('dirDifferenceTensor.png',width = 1000, height = 1000, res= 300)
ggplot(DF, aes(nPC, dirDiff)) + 
  geom_point() + 
  geom_smooth() + 
  theme_bw() + 
  ylab(Average~Difference~from~Symmetric) + 
  xlab(Number~of~PCs)
dev.off()

# for(i in seq_len(ncol(dirMatrix))){
#   for(j in seq_len(ncol(dirMatrix))){
#     if(j != i){
#       if(abs(dirMatrix[i,j]) > abs(dirMatrix[j,i])){
#         dirMatrix[j,i] <- 0
#       } else {
#         dirMatrix[i,j] <- 0
#       }
#     }
#   }
# }

image(dirMatrix, border = FALSE)

X <- t(dirMatrix)
Y <- dirMatrix
Y[20,] <- 0
scTenifoldKnk:::dRegulation(manifoldAlignment(X,Y, 2),1)
