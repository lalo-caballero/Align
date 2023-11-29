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
#' @return A 'list' with:
#' w: Warping function
#' sel: part of y that can be fitted
#' rms: prediction error.
#'
#' @export
#'
#' @examples
#'
#'x <- c(1:10)
#'y <- c(1:10)
#'lambda2 <- 1
#'pow(x, y,lambda2)

pow <- function(x, y, lambda2, lambda1 = 10^6, W = NULL, max_it = 100, min_drms = 1e-6, verbose = FALSE) {
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
    z <- inter$f
    sel <- inter$s
    dg <- inter$g
    z <- z / norm(z, "2")
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
      w <- w + dw
    }
    C <- W %*% G %*% G + lambda2 * Matrix::t(D) %*% D + lambda1 * P
    dw <- as.vector(Matrix::solve(C, (W %*% g * r)))
    w <- w + dw
    diffw <- c(w[2] - w[1], 0.5 * (w[3 : length(w)] - w[1 : (length(w) - 2)]), w[length(w)]-w[(length(w) - 1)])
    diffw <- as.numeric(diffw <= 0)
    P <- methods::as(methods::as(Matrix::Diagonal(x = diffw), "generalMatrix"), "CsparseMatrix")
  }


  if (it == max_it & verbose){
    cat("The computation exceeded the maximum number of iterations")
  }
  return(list(w = w, sel = sel))
}

