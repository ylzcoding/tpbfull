# Log-scale slice update for phi ~ half-Cauchy(0, scale_phi).
# beta_j | nu_j,lambda_j,phi ~ N(0, phi*nu_j/lambda_j).
# This is a full conditional update; no local variables are marginalized.
Gibbs_phi <- function(phi, beta, nu, lambda, scale_phi = 1) {
  p <- length(beta)
  terms <- 2 * log(abs(beta)) + log(lambda) - log(nu) - log(2)
  peak <- max(terms)
  log_quad <- peak + log(sum(exp(terms - peak)))
  log_density <- function(z) {
    (1 - p / 2) * z - exp(log_quad - z) -
      tpb_softplus(2 * (z - log(scale_phi)))
  }

  exp(tpb_slice(log(phi), log_density))
}
