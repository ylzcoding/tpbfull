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

  p <- length(beta_vec)
  prior_type_a <- match.arg(prior_type_a, c("gamma", "hcauchy"))
  prior_type_b <- match.arg(prior_type_b, c("gamma", "hcauchy"))

  r <- 0.5 + b
  s <- 1.5 - a
  z <- (beta_vec^2) / (2 * phi)

  tricomi_U_integral <- function(U_a, U_b, z_val) {
    integrand <- function(t) {
      exp(-z_val * t) * t^(U_a - 1) * (1 + t)^(U_b - U_a - 1)
    }
    integral_result <- tryCatch({
      integrate(integrand, lower = 0, upper = Inf)$value
    }, error = function(e) NA_real_)
    integral_result / gamma(U_a)
  }

  U_vals <- tryCatch({
    gsl::hyperg_U(r, s, z)
  }, error = function(e) {
    return(rep(NaN, p))
  })

  bad_u <- !is.finite(U_vals) | U_vals <= 0
  if (any(bad_u)) {
    U_vals[bad_u] <- vapply(z[bad_u], function(z_val) {
      if (!is.finite(z_val) || z_val < 0) {
        return(NA_real_)
      }
      tricomi_U_integral(r, s, z_val)
    }, numeric(1))
  }

  # if U <= 0，just return -inf to reject this proposal
  if (any(is.nan(U_vals)) || any(U_vals <= 0) || any(is.infinite(U_vals))) {
    return(list(log_posterior = -Inf, log_lik = -Inf))
  }

  term1 <- lgamma(0.5+b) + lgamma(a+b) - lgamma(a) - lgamma(b) - 0.5 * log(2*pi*phi)
  log_lik <- p * term1 + sum(log(U_vals))
  log_prior_a <- if (prior_type_a == "gamma") {
    (s_prior_a - 1) * log(a) - r_prior_a * a
  } else {
    -log(1 + (a / scale_a)^2)
  }
  log_prior_b <- if (prior_type_b == "gamma") {
    (s_prior_b - 1) * log(b) - r_prior_b * b
  } else {
    -log(1 + (b / scale_b)^2)
  }
  # phi ~ Half-Cauchy(0, scale_phi)
  log_prior_phi <- -log(1 + (phi / scale_phi)^2)

  # Unnormalized Log Posterior
  return(list(
    log_posterior = log_lik + log_prior_a + log_prior_b + log_prior_phi,
    log_lik = log_lik
  ))
}
