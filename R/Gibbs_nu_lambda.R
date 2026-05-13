#' Sample nu using the xi augmentation
#'
#' The TPB local scale is represented as
#' beta_j | nu_j, lambda_j, omega_local ~ N(0, omega_local * nu_j / lambda_j),
#' nu_j | xi_j, a follows the inverse-Gaussian update below, and
#' xi_j | nu_j, a ~ Gamma(a, 1 / nu_j).
#'
#' @param p ncol(X)
#' @param phi Local beta variance scale, named phi for compatibility with
#'   Gibbs_beta(). This is omega in the EB sampler and phi * a / b inside
#'   the fully Bayesian sampler.
#' @param beta p*1 coefficient vector
#' @param lambda p*1 local lambda vector
#' @param xi p*1 auxiliary vector
#' @param a TPB shape parameter
#' @return posterior vector nu
#' @export
Gibbs_nu <- function(p, phi, beta, lambda, xi, a) {
  phi <- max(phi, .Machine$double.xmin)
  lambda <- pmax(lambda, .Machine$double.xmin)
  xi <- pmax(xi, .Machine$double.xmin)
  chi <- pmax(beta^2 * lambda / phi, .Machine$double.xmin)
  statmod::rinvgauss(
    p,
    mean = sqrt((chi + 2 * xi) / (2 * a)),
    shape = chi + 2 * xi
  )
}

#' Sample xi in the nu/lambda TPB augmentation
#' @param p ncol(X)
#' @param a TPB shape parameter
#' @param nu p*1 local nu vector
#' @return posterior vector xi
#' @export
Gibbs_xi <- function(p, a, nu) {
  rgamma(p, shape = a, rate = 1 / pmax(nu, .Machine$double.xmin))
}

#' Sample lambda in the nu/lambda TPB augmentation
#' @param p ncol(X)
#' @param b TPB shape parameter
#' @param phi Local beta variance scale, named phi for compatibility with
#'   Gibbs_beta(). This is omega in the EB sampler and phi * a / b inside
#'   the fully Bayesian sampler.
#' @param beta p*1 coefficient vector
#' @param nu p*1 local nu vector
#' @return posterior vector lambda
#' @export
Gibbs_lambda <- function(p, b, phi, beta, nu) {
  phi <- max(phi, .Machine$double.xmin)
  nu <- pmax(nu, .Machine$double.xmin)
  rgamma(
    p,
    shape = b + 0.5,
    rate = b + 0.5 * beta^2 / (phi * nu)
  )
}
