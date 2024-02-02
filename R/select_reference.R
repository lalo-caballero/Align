#' Select reference sample
#'
#' @param data a matrix with all the samples (by row)
#' @param figure 'variance' or 'correlation' figure of merit to select best reference sample
#'
#' @return index of best reference sample
#' @export
#'
#' @examples
#' data <- synthetic_dataset(n_samples = 10,
#'                                  n_peaks = 7,
#'                                  length_out = 100,
#'                                  mov_peaks = 5)
#' select_reference(data)

select_reference <- function(data, figure = 'correlation') {
  data <- as.matrix(data)
  if (figure == 'variance') {
    variance <- apply(data, 1, stats::var)
    idx <- which.max(variance)
  } else if (figure == 'correlation') {
    norms <- sqrt(rowSums(data ^ 2))
    correlation_matrix <- (data %*% t(data)) / (norms %*% t(norms))
    diag(correlation_matrix) <- 0
    mean_correlation <- rowMeans(correlation_matrix, na.rm = TRUE)
    idx <- which.max(mean_correlation)
  } else {
    warning("Not a supported figure of merit, correlation will be used")
    idx <- which.max(mean_correlation)
  }
  return(idx)
}

