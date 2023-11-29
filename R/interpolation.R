#' Linear interpolation.
#'
#' Interpolates a vector.
#'
#' @param t points at which to compute interpolated values.
#' @param y signal to be interpolated.
#'
#' @return a list with:
#' f: the interpolated signal,
#' s: points  of t with  1<= t & t <= length of y
#' g: gradient of y at points t.
#'
#' @export
#' @examples
#' t <- c(1:10)
#' y <- c(1:10)
#' interpolation(t,y)

interpolation<-function(t,y){
  s <- which(1 < t & t < length(y))
  ti <- floor(t[s])
  tr <- t[s] - ti
  g <- y[ti + 1] - y[ti]
  f <- y[ti] + tr * g
  return (list(f = f, s = s, g = g))
}
