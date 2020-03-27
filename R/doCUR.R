doCUR <- function(WT,gKO){
  set.seed(1)
  KO <- CUR(WT, sv = svds(WT, 5))
  C <- KO@C
  C[gKO,] <- 0
  KO <- C %*% KO@U %*% KO@R
  KO <- round(KO, 1)
  return(KO)
}
