#' Compute marginal log posterior for TPB hyperparameters
#' @import gsl
#' @import mvtnorm
#' @param a shape parameter a
#' @param b shape parameter b
#' @param phi global shrinkage parameter
#' @param beta_vec Coefficient vector.
#' @param s_prior_a Prior shape parameter for a.
#' @param r_prior_a Prior rate parameter for a.
#' @param s_prior_b Prior shape parameter for b.
#' @param r_prior_b Prior rate parameter for b.
#' @param prior_type_a prior type for a: "gamma" or "hcauchy"
#' @param prior_type_b prior type for b: "gamma" or "hcauchy"
#' @param scale_a Half-Cauchy scale parameter for a
#' @param scale_b Half-Cauchy scale parameter for b
#' @param scale_phi Half-Cauchy scale parameter for phi.
log_marginal_posterior <- function(a, b, phi, beta_vec,
                                   s_prior_a, r_prior_a,
                                   s_prior_b, r_prior_b,
                                   scale_phi,
                                   prior_type_a = "gamma",
                                   prior_type_b = "gamma",
                                   scale_a = 1,
                                   scale_b = 1) {

  prior_type_a <- match.arg(prior_type_a, c("gamma", "hcauchy"))
  prior_type_b <- match.arg(prior_type_b, c("gamma", "hcauchy"))

  if (!is.finite(a) || !is.finite(b) || !is.finite(phi) ||
      a <= 0 || b <= 0 || phi <= 0) {
    return(list(log_posterior = -Inf, log_lik = -Inf))
  }

  log_lik <- sum(tpb_log_marginal_beta(beta_vec, a, b, phi))
  log_prior_a <- switch(
    prior_type_a,
    gamma = (s_prior_a - 1) * log(a) - r_prior_a * a,
    hcauchy = -tpb_softplus(2 * (log(a) - log(scale_a)))
  )
  log_prior_b <- switch(
    prior_type_b,
    gamma = (s_prior_b - 1) * log(b) - r_prior_b * b,
    hcauchy = -tpb_softplus(2 * (log(b) - log(scale_b)))
  )
  # phi ~ Half-Cauchy(0, scale_phi)
  log_prior_phi <- -tpb_softplus(2 * (log(phi) - log(scale_phi)))

  # Unnormalized Log Posterior
  return(list(
    log_posterior = log_lik + log_prior_a + log_prior_b + log_prior_phi,
    log_lik = log_lik
  ))
}
