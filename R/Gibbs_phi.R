# Log prior density with respect to phi; callers add the log(phi) Jacobian.
tpb_log_phi_prior <- function(log_phi, scale_phi, prior_type_phi = "sqrt_hcauchy") {
  if (identical(prior_type_phi, "hcauchy"))
    return(-tpb_softplus(2 * (log_phi - log(scale_phi))))
  -.5 * log_phi - tpb_softplus(log_phi - 2 * log(scale_phi))
}

# Log-scale slice update for the global variance multiplier.
# beta_j | nu_j,lambda_j,phi ~ N(0, phi*nu_j/lambda_j).
# This is a full conditional update; no local variables are marginalized.
Gibbs_phi <- function(phi, beta, nu, lambda, scale_phi = 1,
                      prior_type_phi = "sqrt_hcauchy") {
  p <- length(beta)
  terms <- 2 * log(abs(beta)) + log(lambda) - log(nu) - log(2)
  peak <- max(terms)
  log_quad <- peak + log(sum(exp(terms - peak)))
  log_density <- function(z) {
    (1 - p / 2) * z - exp(log_quad - z) +
      tpb_log_phi_prior(z, scale_phi, prior_type_phi)
  }

  exp(tpb_slice(log(phi), log_density))
}
