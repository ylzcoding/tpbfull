#' @param p ncol(X)
#' @param a real number
#' @param nu p*1 coefficient vector
Gibbs_xi = function(a, p, nu){return(rgamma(p, shape = a, rate = 1/nu))}


#' @import statmod
#' @param p ncol(X)
#' @param sigmaSq coefficient, real number
#' @param beta p*1 coefficient vector
#' @param lambda p*1 coefficient vector
#' @param xi p*1 coefficient vector, auxiliary variable
#' @export
aug_Gibbs_nu = function(p, phi, beta, lambda, xi){
  chi <- (beta^2 * lambda) / max(phi, 1e-12)
  nu <- statmod::rinvgauss(p, mean = sqrt(chi/2 + xi), shape = chi+2*xi)
  return(nu)
}
