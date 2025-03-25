#' draw posterior lambda
#' @param nu p*1 coefficient vector
#' @param beta p*1 coefficient vector
#' @param omega coefficient, real number
#' @param sigmasq coefficient, real number
#' @param b hyper-parameter, real number
#' @param p ncol(X)
#' @return posterior vector lambda
#' @export
Gibbs_lambda = function(b, p, omega, sigmaSq, beta, nu){
  # lambda_j ~ Gamma(b, b)
  # lambda_j|. ~ Gamma(b+1/2, b+beta_j^2/(2*sigmaSq*omega*nu_j))
  lambda = rgamma(p, shape = b+0.5, rate = b+0.5*beta^2/(sigmaSq*omega*nu))
  return(lambda)
}
