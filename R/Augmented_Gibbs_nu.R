#' @param p ncol(X)
#' @param a hyper-parameter, real number
#' @param nu p*1 coefficient vector
#' @export
Gibbs_xi = function(a, p, nu){return(rgamma(p, shape = a, rate = 1/nu))}


#' @import statmod
#' @param p ncol(X)
#' @param omega coefficient, real number
#' @param sigmasq coefficient, real number
#' @param beta p*1 coefficient vector
#' @param lambda p*1 coefficient vector
#' @param xi p*1 coefficient vector, auxiliary variable
#' @export
aug_Gibbs_nu = function(a, p, omega, sigmaSq, beta, lambda, xi){
  chi = exp(log(beta^2) + log(lambda) - log(sigmaSq) - log(omega))
  nu = statmod::rinvgauss(p, mean = sqrt((chi+2*xi)/(2*a)), shape = chi+2*xi)
  return(nu)
}
