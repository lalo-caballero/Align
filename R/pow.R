#' Penalized Optimized Warping
#'
#' Alignment with POW method.
#'
#' POW like PTW but without a fix polynomial, it will work
#'
#' @param x Signal to fit.
#' @param y Reference signal.
#' @param lambda1 Regularization parameter for second derivative of warp
#' @param lambda2 Penalty for negative values of firsT derivative of warp
#' @param W A 'diagonal matrix' with weights.
#' @param max_it Maximum amount of iterations.
#' @param min_drms minimum difference in absolute value of RMS.
#' @param verbose by default 'FALSE'to privide information in the console
#'
#' @return warp
#'
#' @export
#'
#' @examples
#'
#'x <- c(1:10)
#'y <- c(1:10)
#'lambda2 <- 1
#'pow(x, lambda2, y)

pow <- function(x, lambda2, y, lambda1 = 10^6, W = NULL, max_it = 100, min_drms = 1e-6, verbose = FALSE) {
  require(Matrix)
  m<-max(length(x),length(y))
  if (is.null(W)){
    W <- methods::as(methods::as(Matrix::Diagonal(m, 1), "generalMatrix"), "CsparseMatrix")
  } else {
    W <- methods::as(methods::as(W, "generalMatrix"), "CsparseMatrix")
  }
  t <- 1:m
  w <- 1:m
  rms_old <- 0
  x <- x / norm(x, "2")
  y <- y / norm(y, "2")
  I <- Matrix::diag(length(y))
  D <- diff(I, differences = 2)
  D <- methods::as(D, "dgCMatrix")
  for (it in 1 : max_it) {
    r <- rep(0, m)
    g <- rep(0, m)
    inter <- interpolation(w, x)
    z <- inter$f / norm(inter$f, "2")
    sel <- inter$s
    dg <- inter$g
    g[sel] <- dg
    G <- methods::as(methods::as(Matrix::Diagonal(x = g), "generalMatrix"), "CsparseMatrix")
    r[sel] <- y[sel] - z
    rms <- sqrt(sum(r ^ 2) / m)
    drms <- abs((rms - rms_old) / (rms + 1e-10))
    if (drms < min_drms) {
      break
    }
    rms_old <- rms
    if (it == 1) {
      C <- W %*% G %*% G + lambda2 * Matrix::t(D) %*% D
      dw <- as.vector(Matrix::solve(C, (W %*% g * r)))
      w <- w + dw
      diffw <- c(w[2] - w[1], 0.5 * (w[3 : length(w)] - w[1 : (length(w) - 2)]), w[length(w)]-w[(length(w) - 1)])
      diffw <- as.numeric(diffw <= 0)
      P <- methods::as(methods::as(Matrix::Diagonal(x = diffw), "generalMatrix"), "CsparseMatrix")
      w <- w - dw
    }
    C <- W %*% G %*% G + lambda2 * Matrix::t(D) %*% D + lambda1 * P
    dw <- as.vector(Matrix::solve(C, (W %*% g * r), tol = 1e-20))
    w <- w + dw
    diffw <- c(w[2] - w[1], 0.5 * (w[3 : length(w)] - w[1 : (length(w) - 2)]), w[length(w)]-w[(length(w) - 1)])
    diffw <- as.numeric(diffw <= 0)
    P <- methods::as(methods::as(Matrix::Diagonal(x = diffw), "generalMatrix"), "CsparseMatrix")
  }


  if (it == max_it & verbose){
    cat("The computation exceeded the maximum number of iterations")
  }
  return(w)
}


#' Apply POW
#'
#' @param X A matrix with the samples to align
#' @param y The reference signal
#' @param lambdas the lambda to be used for the alignment of each sample
#' @param max_it maximum iterations to be done in the POW alignment
#'
#' @return a list with:
#' warped_samples
#' warps
#' @export
#'
#' @examples
apply_pow <- function(X, lambdas, y, max_it = 1000){

  n_samples <- nrow(X)
  XW <- X * 0
  for (i in 1:n_samples){
    w <- pow(X[i, ], lambdas[i], y, max_it = max_it)
    interp <- interpolation(w, X[i, ])
    xw <- interp$f
    sel <- interp$s
    XW[i,sel] <- xw
    XW[i,-sel] <- NA
  }
  return (XW)
}
