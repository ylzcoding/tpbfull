#' @param phi global shrinkage parameter
#' @param beta p*1 coefficient vector
#' @param lambda p*1 coefficient vector
#' @param a shape parameter from the prior nu ~ Gamma(a, 1)
#' @return p*1 vector nu_new
#' @importFrom GIGrvg rgig
#' @export
Gibbs_nu <- function(phi, beta, lambda, a) {
  p <- length(beta)
  nu = c()
  for (d in 1:p){
    chi = max(exp(log(beta[d]^2) + log(lambda[d]) - log(phi)),
              10^(-306)) # 307 is smallest integer that didn't result in the code breaking here
    nu = c(nu, GIGrvg::rgig(1, lambda = a-0.5, chi = chi, psi = 2))

  }
  return(nu)
}
