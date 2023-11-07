#' Penalized Optimized Warping
#'
#' Alignment with POW method.
#'
#' POW like PTW but without a fix polynomial, it will work
#'
#' @param x Signal to align.
#' @param y A reference "signal".
#' @param lambda1 Regularization parameter for second derivative of warp
#' @param lambda2 Penalty for negative values of firs derivative (warp)
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
#'pow(x, y)

pow <- function(x, y, lambda1, lambda2 = 10^6, W=NULL, max_it=1000, min_drms=1e-5,verbose=FALSE){
  #library(Matrix)
  m<-max(length(x),length(y))
  if (is.null(W)){
    W<-as(as(Matrix::Diagonal(m,1),"generalMatrix"),"CsparseMatrix")
  } else {
    W<-as(as(W,"generalMatrix"),"CsparseMatrix")
  }
  t<-1:m
  w<-1:m
  rms<-rep(NA,1000)
  #rms[1]<-0
  #drms<-rep(NA,1000)
  x<-x/norm(x,"2")
  y<-y/norm(y,"2")
  I<-Matrix::diag(length(y))
  D<-diff(I)
  D<-as(D,"dgCMatrix")
  for (it in 1:max_it){
    r<-rep(0,m)
    g<-rep(0,m)
    inter<-interpolation(w,x)
    z<-inter$f
    sel<-inter$s
    dg<-inter$g
    if (length(z)<20){
      break
    }
    z<-z/norm(z,"2")
    g[sel]<-dg
    G<-as(as(Matrix::Diagonal(x=g),"generalMatrix"),"CsparseMatrix")
    r[sel]<-y[sel]-z
    rms[it] <- sqrt(sum(r^2) / m)
    #drms[it]<-abs((rms[it]-rms[it-1])/(rms[it]+1e-10))
    # if (drms[it] < min_drms){
    #   break
    # }
    dw<-Matrix::solve((W %*% G %*% G + lambda1 * Matrix::t(D) %*% D),(W %*% g * r))
    w<-w+as.vector(dw)
  }
  #derivative of w
  #dw*y
  if (it== max_it & verbose){
    cat("The computation exceeded the maximum number of iterations")
  }
  return(list(w=w,rms=rms))
}

