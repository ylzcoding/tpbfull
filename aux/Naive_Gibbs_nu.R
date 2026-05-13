#' @import GIGrvg
#' @param p ncol(X)
#' @param a numeric, one of the hyperparameters
#' @param phi coefficient, real number
#' @param sigmasq coefficient, real number
#' @param beta p*1 coefficient vector
#' @param lambda p*1 coefficient vector
#' @export
Naive_Gibbs_nu = function(a, p, omega, beta, lambda){
  # nu_j|a,1 ~ Gamma(a, 1)
  nu = c()
  for (d in 1:p){
    chi = max(exp(log(beta[d]^2) + log(lambda[d]) - log(omega)),
               10^(-306)) # 307 is smallest integer that didn't result in the code breaking here
    nu = c(nu, GIGrvg::rgig(1, lambda = a-0.5, chi = chi, psi = 2*a))

  }
  return(nu)
}
