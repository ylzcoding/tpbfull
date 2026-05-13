#' @import statmod
#' @param p ncol(X)
#' @param phi coefficient, real number
#' @param sigmasq coefficient, real number
#' @param beta p*1 coefficient vector
#' @param lambda p*1 coefficient vector
#' @param xi p*1 coefficient vector, auxiliary variable
#' @export
aug_Gibbs_nu = function(p, omega, beta, lambda, xi, a){
  chi = exp(log(beta^2) + log(lambda) - log(omega))
  nu = statmod::rinvgauss(p, mean = sqrt((chi+2*xi)/(2*a)), shape = chi+2*xi)
  return(nu)
}
