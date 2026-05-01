#' Sample error variance sigmaSq
#'
#' @param n Number of observations.
#' @param X n*p design matrix.
#' @param y n*1 response vector.
#' @param beta p*1 coefficient vector.
#' @return A scalar draw of sigmaSq.
Gibbs_sigmaSq = function(n, X, y, beta){
  return(1/rgamma(1, shape = n/2+3/2, rate = 1/2 + 0.5*(t(y - X%*%beta)%*%(y - X%*%beta))))
}
