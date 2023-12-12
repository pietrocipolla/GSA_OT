bw_OT_indices <- function(x, y, M) {
  # Transform data to match the respective formats
  y <- y %>%
    pivot_wider(names_from = t) %>%
    select(where(~ any(. != 0)))
  y_cols <- colnames(y)
  x <- x %>%
    mutate(run = as.factor(row_number()))
  
  # Remove all the rows that did not reach the end of the simulation
  xy <- x %>%
    inner_join(y, by = "run") %>%
    filter(!if_any(everything(), is.na))
  
  # Retrieve the cleaned data
  x <- xy %>%
    select(-all_of(y_cols))
  y <- xy %>%
    select(all_of(y_cols)) %>%
    select(-run) %>%
    as.matrix()
  
  # Compute the statistics for the unconditioned distribution
  my <- colMeans(y)
  Cy <- cov(y)
  traceCy <- sum(diag(Cy))
  Ry <- sqrtm(Cy)
  
  # Retrieve values useful for the algorithm
  N <- dim(x)[1]
  K <- dim(x)[2]
  
  # Define the partitions for the estimation
  partitions_indices <- floor(seq(from = 1, to = N + 1, length.out = M + 1))
  
  # Initialize the result matrices
  W <- matrix(nrow = M, ncol = K)
  Adv <- matrix(nrow = M, ncol = K)
  Diff <- matrix(nrow = M, ncol = K)
  
  # Evaluate the upper bound of the indices
  V <- 2 * sum(diag(Cy))
  
  # Run the algorithm for each variable
  for (k in seq_len(K)) {
    # Get the ranking of the current variable samples
    ord <- rank(x[, k], ties.method = "random")
    
    ret <- bind_rows(map(seq_len(M), bw, 
                         partitions_indices, y, ord,
                         my, Cy, traceCy, Ry))
    
    Adv[, k] <- ret$Adv
    Diff[, k] <- ret$Diff
    W[, k] <- ret$W
  }
  
  indices <- matrix(nrow = 3, ncol = K)
  
  indices[1, ] <- colMeans(W) / V
  indices[2, ] <- colMeans(Adv) / V
  indices[3, ] <- colMeans(Diff) / V
  colnames(indices) <- colnames(x)
  rownames(indices) <- c("BW", "Adv", "Diff")
  
  return(indices)
}
