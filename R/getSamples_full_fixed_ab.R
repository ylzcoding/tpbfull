#' Fully Bayesian Gibbs Sampler for the TPB Normal Model with Fixed a and b
#' Uses sigma^2 ~ InvGamma(1.5, 0.5) and sqrt(phi) ~ half-Cauchy(0, scale_phi),
#' matching the main TPB sampler.
#' @import GIGrvg
#' @param num_burnin The number of MCMC iterations to discard as burn-in.
#' @param num_samples The number of MCMC iterations to save after the burn-in period.
#' @param fix_phi Logical. If TRUE, phi is fixed to the provided phi_val without sampling.
#' @param phi_val Numeric. The fixed value for phi (if fix_phi = TRUE) or the initial starting value.
#' @param scale_phi Half-Cauchy scale of the quantity selected by prior_type_phi.
#' @param prior_type_phi "sqrt_hcauchy" for a prior on sqrt(phi), or "hcauchy"
#'   for the earlier prior on phi itself.
#' @param m0 Optional prior mean effective signal count; calibrates scale_phi.
#' @param calibration_points Number of Sobol points used for calibration.
#' @param scale_with_sigma If TRUE, the coefficient variance is sigmaSq*phi*nu/lambda.
#'   This supports the horseshoe PV prior: set a = b = 0.5, prior_type_phi =
#'   "sqrt_hcauchy", and scale_phi = m0/(p-m0)/sqrt(n) for standardized columns.
#'   Leave m0 = NULL for this PV reference rule; m0 instead invokes the full-prior
#'   mean calibration. Returned phi is the multiplier before sigmaSq.
#' @param store_locals Save nu and lambda draws. FALSE saves only beta and scalars.
#' @return list of posterior samples
#' @export

tpb_fullyBayes_fixab <- function(X, y, a, b, num_burnin, num_samples,
                                 woodbury = TRUE, diagX = FALSE,
                                 fix_phi = FALSE, phi_val = NULL,
                                 scale_phi = 1,
                                 prior_type_phi = c("sqrt_hcauchy", "hcauchy"),
                                 m0 = NULL, calibration_points = 262144,
                                 scale_with_sigma = FALSE, store_locals = TRUE) {
  prior_type_phi <- match.arg(prior_type_phi)
  if (fix_phi && is.null(phi_val)) stop("Supply phi_val when fix_phi = TRUE.")
  if (fix_phi && !is.null(m0)) stop("m0 cannot calibrate a fixed phi.")
  calibration <- NULL
  if (!is.null(m0)) {
    calibration <- tpb_calibrate_phi(X, m0,
      list(a = a, b = b, prior_type_phi = prior_type_phi),
      sigmaSq = if (scale_with_sigma) 1 else NULL,
      n_points = calibration_points, model = "fixed")
    scale_phi <- calibration$scale_phi
  }
  n <- nrow(X)
  p <- ncol(X)
  num_total <- num_burnin + num_samples
  if (is.null(scale_phi) || !is.finite(scale_phi) || scale_phi <= 0) {
    scale_phi <- 1
  }
  
  beta_samples <- matrix(NA, nrow = num_samples, ncol = p)
  nu_samples <- lambda_samples <-
    if (store_locals) matrix(NA, nrow = num_samples, ncol = p) else NULL
  sigmaSq_samples <- numeric(num_samples)
  phi_samples <- numeric(num_samples)
  
  ### Initial values
  sigmaSq_cur <- var(y)
  phi_cur <- if (!is.null(phi_val)) phi_val else
    if (prior_type_phi == "sqrt_hcauchy") scale_phi^2 else scale_phi
  nu_cur <- rgamma(p, 1, 1)
  lambda_cur <- rgamma(p, 1, 1)
  
  idx_store <- 0
  
  for (i in 1:num_total) {
    variance_scale <- phi_cur * if (scale_with_sigma) sigmaSq_cur else 1
    beta_cur <- Gibbs_beta(
      X = X, y = y, sigmaSq = sigmaSq_cur,
      phi = variance_scale, nu = nu_cur, lambda = lambda_cur,
      woodbury = woodbury, diagX = diagX
    )

    for (d in 1:p){
      chi <- max(exp(2 * log(abs(beta_cur[d])) + log(lambda_cur[d]) - log(variance_scale)),
                 10^(-306)) # 307 is smallest integer that didn't result in the code breaking here
      nu_cur[d] <- GIGrvg::rgig(1, lambda = a-0.5, chi = chi, psi = 2)
    }
    rate_lam <- 1 + beta_cur^2 / (2*variance_scale*nu_cur)
    lambda_cur <- rgamma(p, shape = b+0.5, rate  = rate_lam)
    
    if (scale_with_sigma) {
      residual <- y - drop(X %*% beta_cur)
      prior_quad <- sum(beta_cur^2 * lambda_cur / (phi_cur * nu_cur))
      sigmaSq_cur <- 1 / rgamma(1, (n+p+3)/2,
                                .5 + .5*(sum(residual^2) + prior_quad))
    } else {
      sigmaSq_cur <- Gibbs_sigmaSq(n, X, y, beta_cur)
    }
    if (!fix_phi) {
      beta_scaled <- if (scale_with_sigma) beta_cur / sqrt(sigmaSq_cur) else beta_cur
      phi_cur <- Gibbs_phi(phi_cur, beta_scaled, nu_cur, lambda_cur,
                           scale_phi, prior_type_phi)
    }
    
    # store after burn-in
    if(i > num_burnin) {
      idx_store <- idx_store + 1
      beta_samples[idx_store, ] <- beta_cur
      if (store_locals) {
        lambda_samples[idx_store, ] <- lambda_cur
        nu_samples[idx_store, ] <- nu_cur
      }
      sigmaSq_samples[idx_store] <- sigmaSq_cur
      phi_samples[idx_store] <- phi_cur
    }
  }
  return(
    list(beta = beta_samples,
         lambda = lambda_samples,
         nu = nu_samples,
         sigmaSq = sigmaSq_samples,
         phi = phi_samples, scale_phi = scale_phi, prior_type_phi = prior_type_phi,
         scale_with_sigma = scale_with_sigma, calibration = calibration)
  )
}
