#' @param p ncol(X)
#' @param w numeric, auxiliary variable for sampling from half-cauchy dist.
#' @param lambda p*1 coefficient vector
#' @param nu p*1 coefficient vector
#' @param beta p*1 coefficient vector
#' @param sigmaSq coefficient, real number
#' @export
Gibbs_omega = function(p, nu, lambda, w, sigmaSq, beta){
  return(1/rgamma(1, shape = (p+1)/2, rate = 1/w+1/(2*sigmaSq)*sum(beta^2*lambda/nu)))
}


#' @param omega_hyper numeric, omega ~ half-cauchy(0, omega_hyper), default value: 1
#' @export
Gibbs_w = function(omega, omega_hyper = 1){return(1/rgamma(1, 1, 1/omega+1/omega_hyper))}

#' @param p ncol(X)
#' @param n nrow(X)
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param nu p*1 coefficient vector
#' @param beta p*1 coefficient vector
#' @param omega coefficient, real number
#' @param lambda p*1 coefficient vector
#' @export
Gibbs_sigmaSq = function(n, p, y, X, beta, omega, nu, lambda){
  eta = omega*nu/lambda
  return(1/rgamma(1, shape = (n+p)/2, rate = 0.5*(t(y - X%*%beta)%*%(y - X%*%beta)+t(beta)%*%diag(1/eta)%*%beta)))
}
