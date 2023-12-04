#' Select reference sample
#'
#' @param data a matrix with all the samples (by row)
#' @param figure 'variance' or 'correlation' figure of merit to select best reference sample
#'
#' @return index of best reference sample
#' @export
#'
#' @examples

select_reference <- function(data, figure = 'variance') {
  data <- as.matrix(data)
  n_samples <- nrow(data)
  norms <- sqrt(rowSums(data^2))
  correlation_matrix <- (data %*% t(data)) / (norms %*% t(norms))
  diag(correlation_matrix) <- 0
  mean_correlation <- rowMeans(correlation_matrix, na.rm = TRUE)
  variance <- apply(data, 1, stats::var)
  if (figure == 'variance') {
    idx <- which.max(variance)
  } else if (figure == 'correlation') {
    idx <- which.max(mean_correlation)
  } else {
    warning("Not a supported figure of merit, variance will be used")
    idx <- which.max(variance)
  }
  return(idx)
}

