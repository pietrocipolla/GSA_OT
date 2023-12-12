ot_indices_1d <- function(x, y, c = "L2", M = 25) {
  # Retrieve x ranks
  ord <- rank(x)
  
  # Sort y
  y_sort <- sort(y)
  
  # Get the number of realizations
  N <- length(x)
  
  # Build the partitions. Each partition has ~ the same number of elements
  partitions_indices <- floor(seq(from = 1, to = N + 1, length.out = M + 1))
  partitions <- list()
  
  for (m in seq_len(M)) {
    partitions[[m]] <- partitions_indices[m]:(partitions_indices[m+1] - 1)
    partitions[[m]] <- which(ord %in% partitions[[m]])
  }
  
  # Get the number of elements in each partition
  n <- diff(partitions_indices)
  
  W <- future_map_dbl(partitions, optimal_transport_1d, y, y_sort,
                      .options = furrr_options(seed = 777))
  
  return((W %*% n) / (2 * var(y) * N))
}
