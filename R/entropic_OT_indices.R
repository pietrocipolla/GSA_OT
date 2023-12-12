entropic_OT_indices <- function(x, y, M, metric = "L2") {
  C = as.matrix(dist(y, method = "euclidean"))
  N <- dim(x)[1]
  K <- dim(x)[2]
  
  partitions_indices <- floor(seq(from = 1, to = N + 1, length.out = M + 1))
  
  W <- matrix(nrow = M, ncol = K)
  
  V <- 2 * sum(diag(cov(y)))
  
  for (k in seq_len(K)) {
    ord <- rank(x[, k])
    
    for (m in seq_len(M)) {
      partition <- partitions_indices[m]:(partitions_indices[m + 1] - 1)
      ii <- which(ord %in% partition)
      
      out <- sinkhorn(C[ii, ], 10000, 0.01)
      W[m, k] <- out$W22
    }
  }
  
  return(W / V)
}