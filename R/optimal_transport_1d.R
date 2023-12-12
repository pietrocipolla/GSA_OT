# y and y_sort are both included for computational reasons

optimal_transport_1d <- function(partition, y, y_sort) {
  # Extract the conditioned empirical CDF quantiles
  yy <- sort(y[partition])
  
  # Set the scaling parameters
  N <- length(y)
  n <- length(yy)
  
  # Expand the empirical CDF quantiles to match with y length
  yy <- yy[floor(seq(1/n, 1, length.out = N) * n + 0.5)]
  
  W <- mean((y_sort - yy)^2)
}
