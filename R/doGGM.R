doGGM <- function(WT, gKO){
  diag(WT) <- 1
  gKO <- which(rownames(WT) %in% gKO)
  n <- nrow(WT)
  if (index_rm == 1){
    per_net_mat <- WT
  } else {
    pet_net_mat <- WT[-gKO,-gKO]
  }
  
  ## get inverse of submatrix
  q11 <- per_net_mat[1, 1]
  q1 <- per_net_mat[1, -1]
  q1t <- per_net_mat[-1, 1]
  net_mat_inv <- per_net_mat[-1,-1] - tcrossprod(q1t, q1)/q11
  
  ## add matrix names
  if (is.null(rownames(WT))){
    return(net_mat_inv)
  } else{
    name_mat <- rownames(WT)
    rownames(net_mat_inv) <- name_mat[-gKO]
    colnames(net_mat_inv) <- name_mat[-gKO]
    return(net_mat_inv)
  }
}