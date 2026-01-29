#' draw posterior lambda
#' @param nu p*1 coefficient vector
#' @param beta p*1 coefficient vector
#' @param omega coefficient, real number
#' @param b hyper-parameter, real number
#' @param p ncol(X)
#' @return posterior vector lambda

Gibbs_lambda = function(b, p, phi, beta, nu){
  # lambda_j|b ~ Gamma(b, 1)
  lambda = rgamma(p, shape = b+0.5, rate = 1+0.5*beta^2/(phi*nu))
  return(lambda)
}
