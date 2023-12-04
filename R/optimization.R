#' Compute validation error.
#'
#' Computes validation error and checks for time inversion in warp function.
#'
#' @param X Matrix with data to be aligned (one sample per row).
#' @param y Reference sample.
#' @param W A matrix with of weights.
#' @param fom figure of merit to use ('rms' or 'cor').
#' @param iv indexes of the points used to validate.
#' @param lambdas a vector with the penalties to test the POW.
#'
#' @return a list with 2 values 'e_ix' and 'ti_ix'.
#' @export
#'
#' @examples

compute_val_error_pow <- function(X, y, W, fom, iv, lambdas) {
  ti_ix <- matrix(0, nrow = nrow(X), ncol = length(lambdas))
  e_ix <- matrix(0, nrow = nrow(X), ncol = length(lambdas))
  s <- 0

  for (samples in 1:nrow(X)) {
    s <- s + 1
    x <- X[, s]
    e <- rep(0,length(lambdas))
    ti <- rep(0,length(lambdas))
    k2 <- 0

    for (lambda2 in lambdas) {
      k2 <- k2 + 1
      w <- pow(x, y, lambda2, W = W)$w
      diff_w <- diff(w)
      ti[k2] <- any(diff_w < 0)
      # compute error
      result <- interpolation(w, x)
      xw <- result$f
      sel <- result$s
      sel_x <- intersect(sel, iv)
      sel_xw <- match(sel_x, sel)

      y_norm <- y[sel_x] / norm(y[sel_x],"2")
      xw_norm <- xw[sel_xw] / norm(xw[sel_xw],"2")

      if (fom == 'rms') {
        e[k2] <- sqrt(sum((y_norm - xw_norm) * (y_norm - xw_norm)) / length(xw_norm))
      } else if (fom == 'cor') {
        e[k2] <- 1 - sum(y_norm * xw_norm)
      }
    }
    # accumulate time inversion check and error
    ti_ix[s, ] <- ti
    e_ix[s, ] <- e
  }

  return(list(e_ix = e_ix, ti_ix = ti_ix))
}


#' Optimize parameters POW
#'
#' @param n_samples number of samples or signals to be used
#' @param params parametros
#' @param ti_ix matrix with time inversions
#' @param e_ix error matrix
#'
#' @return a list with two values best_params and params_ix
#' @export
#'
#' @examples
optimize_params_pow <- function(n_samples, params, ti_ix, e_ix) {
  params_ix <- matrix(rep(params, each = n_samples), nrow = n_samples, byrow = TRUE)
  e_ix[ti_ix == 1] <- NA
  params_ix[ti_ix == 1] <- NA
  pos_e_ix <- apply(e_ix, 1, which.min)
  best_params <- params[pos_e_ix]

  return(list(best_params = best_params, params_ix = params_ix))
}

