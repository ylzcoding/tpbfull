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
#' @param prior_type_a prior type for a: "gamma", "hcauchy", or "uniform"
#' @param prior_type_b prior type for b: "gamma", "hcauchy", or "uniform"
#' @param scale_a Half-Cauchy scale parameter for a
#' @param scale_b Half-Cauchy scale parameter for b
#' @param lower_a,upper_a Uniform prior support for a
#' @param lower_b,upper_b Uniform prior support for b
#' @param scale_phi Half-Cauchy scale parameter for phi.
log_marginal_posterior <- function(a, b, phi, beta_vec,
                                   s_prior_a, r_prior_a,
                                   s_prior_b, r_prior_b,
                                   scale_phi,
                                   prior_type_a = "gamma",
                                   prior_type_b = "gamma",
                                   scale_a = 1,
                                   scale_b = 1,
                                   lower_a = 0.01,
                                   upper_a = 10,
                                   lower_b = 0.01,
                                   upper_b = 10) {

  p <- length(beta_vec)
  prior_type_a <- match.arg(prior_type_a, c("gamma", "hcauchy", "uniform"))
  prior_type_b <- match.arg(prior_type_b, c("gamma", "hcauchy", "uniform"))

  if (!is.finite(a) || !is.finite(b) || !is.finite(phi) ||
      a <= 0 || b <= 0 || phi <= 0) {
    return(list(log_posterior = -Inf, log_lik = -Inf))
  }
  if ((prior_type_a == "uniform" &&
       (!is.finite(lower_a) || !is.finite(upper_a) || lower_a <= 0 || lower_a >= upper_a)) ||
      (prior_type_b == "uniform" &&
       (!is.finite(lower_b) || !is.finite(upper_b) || lower_b <= 0 || lower_b >= upper_b))) {
    return(list(log_posterior = -Inf, log_lik = -Inf))
  }

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

  bad_u <- is.na(U_vals) | !is.finite(U_vals) | U_vals <= 0
  if (any(bad_u)) {
    U_vals[bad_u] <- vapply(z[bad_u], function(z_val) {
      if (!is.finite(z_val) || z_val < 0) {
        return(NA_real_)
      }
      tricomi_U_integral(r, s, z_val)
    }, numeric(1))
  }

  # if U <= 0，just return -inf to reject this proposal
  if (any(is.na(U_vals)) || any(!is.finite(U_vals)) || any(U_vals <= 0)) {
    return(list(log_posterior = -Inf, log_lik = -Inf))
  }

  term1 <- lgamma(0.5+b) + lgamma(a+b) - lgamma(a) - lgamma(b) - 0.5 * log(2*pi*phi)
  log_lik <- p * term1 + sum(log(U_vals))
  log_prior_a <- switch(
    prior_type_a,
    gamma = (s_prior_a - 1) * log(a) - r_prior_a * a,
    hcauchy = -log(1 + (a / scale_a)^2),
    uniform = if (a >= lower_a && a <= upper_a) {
      -log(upper_a - lower_a)
    } else {
      -Inf
    }
  )
  log_prior_b <- switch(
    prior_type_b,
    gamma = (s_prior_b - 1) * log(b) - r_prior_b * b,
    hcauchy = -log(1 + (b / scale_b)^2),
    uniform = if (b >= lower_b && b <= upper_b) {
      -log(upper_b - lower_b)
    } else {
      -Inf
    }
  )
  # phi ~ Half-Cauchy(0, scale_phi)
  log_prior_phi <- -log(1 + (phi / scale_phi)^2)

  # Unnormalized Log Posterior
  return(list(
    log_posterior = log_lik + log_prior_a + log_prior_b + log_prior_phi,
    log_lik = log_lik
  ))
}
