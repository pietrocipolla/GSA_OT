bw <- function(m, partitions_indices, y, ord, my, Cy, traceCy, Ry) {
  # Get the current partition indices for y
  partition <- partitions_indices[m]:(partitions_indices[m + 1] - 1)
  ii <- which(ord %in% partition)
  
  yc <- y[ii,]
  
  # Get the conditioned statistics
  mc <- colMeans(yc)
  Cc <- cov(yc)
  
  Adv <- sum((my - mc)^2)
  Diff <- traceCy + sum(diag(Cc)) - 2 * tracesqrtm(Ry %*% Cc %*% Ry)
  
  W <- Adv + Diff
  
  return(data.frame(m, Adv, Diff, W))
}
