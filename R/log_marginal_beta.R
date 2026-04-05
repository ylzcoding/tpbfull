#' compute logp(beta_i|a,b,phi) (integrate out nu and lambda)
#' @import gsl
#' @import mvtnorm
#' @param a shape parameter a
#' @param b shape parameter b
#' @param phi global shrinkage parameter
#' @param beta_vec
#' @param s_prior_a, prior shape parameter for a
#' @param r_prior_a, prior rate parameter for a
#' @param s_prior_b, prior shape parameter for b
#' @param r_prior_b, prior rate parameter for b
#' @param scale_phi
log_marginal_posterior <- function(a, b, phi, beta_vec,
                                   s_prior_a, r_prior_a,
                                   s_prior_b, r_prior_b,
                                   scale_phi) {

  p <- length(beta_vec)

  r <- 0.5 + b
  s <- 1.5 - a
  z <- (beta_vec^2) / (2 * phi)

  U_vals <- tryCatch({
    gsl::hyperg_U(r, s, z)
  }, error = function(e) {
    return(rep(NaN, p))
  })

  # if U <= 0，just return -inf to reject this proposal
  if (any(is.nan(U_vals)) || any(U_vals <= 0) || any(is.infinite(U_vals))) {
    return(-Inf)
  }

  term1 <- lgamma(0.5+b) + lgamma(a+b) - lgamma(a) - lgamma(b) - 0.5 * log(2*pi*phi)
  log_lik <- p * term1 + sum(log(U_vals))
  log_prior_a <- (s_prior_a - 1) * log(a) - r_prior_a * a
  log_prior_b <- (s_prior_b - 1) * log(b) - r_prior_b * b
  # phi ~ Half-Cauchy(0, scale_phi)
  log_prior_phi <- -log(1 + (phi / scale_phi)^2)

  # Unnormalized Log Posterior
  return(log_lik + log_prior_a + log_prior_b + log_prior_phi)
}
