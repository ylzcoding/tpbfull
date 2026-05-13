#' @param p ncol(X)
#' @param w numeric, auxiliary variable for sampling from half-cauchy dist.
#' @param lambda p*1 coefficient vector
#' @param nu p*1 coefficient vector
#' @param beta p*1 coefficient vector
#' @param phi coefficient, real number
#' @export
Gibbs_phi = function(p, nu, lambda, w, beta, prshapeo, prrateo){
  return(1/rgamma(1, shape = p/2 + 1 + prshapeo, rate = sum(beta^2*lambda/nu)/2 + prrateo))
}


#' @param phi_hyper numeric, phi ~ half-cauchy(0, phi_hyper), default value: 1
#' @export
Gibbs_w = function(omega, omega_hyper){return(1/rgamma(1, 1, 1/omega+1/omega_hyper))}

#' @param p ncol(X)
#' @param n nrow(X)
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param nu p*1 coefficient vector
#' @param beta p*1 coefficient vector
#' @param phi coefficient, real number
#' @param lambda p*1 coefficient vector
#' @export
Gibbs_sigmaSq = function(n, p, y, X, beta, prshapes, prrates){
  return(1/rgamma(1, shape = n/2 + 1 + prshapes, rate = sum((y - X%*%beta)^2)/2 + prrates))
}
