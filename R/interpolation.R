#' Interpolation.
#'
#' Interpolates a vector.
#'
#' @param w points at which to compute interpolated values.
#' @param x signal to be interpolated.
#' @param return if TRUE it will return f,s, and g; if FALSE it will return the warped x
#'
#' @return a list with:
#' f: the interpolated signal,
#' s: points  of t with  1<= t & t <= length of y
#' g: gradient of y at points t.
#' or the warped x
#'
#'
#' @export
#' @examples
#' w <- c(1:10)
#' x <- c(1:10)
#' interpolation(w,x)

interpolation <- function (w, x, return = TRUE){
  s <- which(1 < w & w < length(x))
  ti <- floor(w[s])
  tr <- w[s] - ti
  g <- x[ti + 1] - x[ti]
  f <- x[ti] + tr * g
  x[s] <- f
  x[-s] <- NA
  if (return){
    return (list(f = f, s = s, g = g))
  } else{
    return(x)
  }

}
