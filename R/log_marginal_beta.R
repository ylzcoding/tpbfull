#' Compute marginal log posterior for TPB hyperparameters
#' @import gsl
#' @import mvtnorm
#' @param a shape parameter a
#' @param b shape parameter b
#' @param phi global shrinkage parameter
#' @param beta_vec Coefficient vector.
#' @param log_abs_beta Optional log absolute coefficients, preserving tiny values.
#' @param hyper_params Hyperprior specification used by fullGibbs().
log_marginal_posterior <- function(a, b, phi, beta_vec, hyper_params,
                                   log_abs_beta = log(abs(beta_vec))) {
  if (!is.finite(a) || !is.finite(b) || !is.finite(phi) ||
      a <= 0 || b <= 0 || phi <= 0) {
    return(list(log_posterior = -Inf, log_lik = -Inf))
  }

  log_lik <- sum(tpb_log_marginal_beta(beta_vec, a, b, phi, log_abs_beta))
  log_prior <- tpb_log_hyperprior(log(c(a, b, phi)), hyper_params)

  # Unnormalized Log Posterior
  return(list(
    log_posterior = log_lik + log_prior,
    log_lik = log_lik
  ))
}
