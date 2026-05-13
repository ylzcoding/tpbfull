#' draw posterior lambda
#' @param nu p*1 coefficient vector
#' @param beta p*1 coefficient vector
#' @param phi coefficient, real number
#' @param sigmasq coefficient, real number
#' @param b hyper-parameter, real number
#' @param p ncol(X)
#' @return posterior vector lambda
#' @export
Gibbs_lambda = function(b, p, omega, beta, nu){
  # lambda_j|b ~ Gamma(b, 1)
  lambda = rgamma(p, shape = b+0.5, rate = b+0.5*beta^2/(omega*nu))
  return(lambda)
}
