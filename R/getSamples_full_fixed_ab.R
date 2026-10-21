#' Fully Bayesian Gibbs Sampler for the TPB Normal Model with Fixed a and b
#' Uses sigma^2 ~ InvGamma(1.5, 0.5) and phi ~ half-Cauchy(0, scale_phi),
#' matching the main TPB sampler.
#' @import GIGrvg
#' @param num_burnin The number of MCMC iterations to discard as burn-in.
#' @param num_samples The number of MCMC iterations to save after the burn-in period.
#' @param fix_phi Logical. If TRUE, phi is fixed to the provided phi_val without sampling.
#' @param phi_val Numeric. The fixed value for phi (if fix_phi = TRUE) or the initial starting value.
#' @param scale_phi Half-Cauchy scale for phi (not sqrt(phi)).
#' @return list of posterior samples
#' @export

tpb_fullyBayes_fixab <- function(X, y, a, b, num_burnin, num_samples,
                                 woodbury = TRUE, diagX = FALSE,
                                 fix_phi = FALSE, phi_val = NULL,
                                 scale_phi = 1) {
  n <- nrow(X)
  p <- ncol(X)
  num_total <- num_burnin + num_samples
  if (is.null(scale_phi) || !is.finite(scale_phi) || scale_phi <= 0) {
    scale_phi <- 1
  }
  
  beta_samples <- matrix(NA, nrow = num_samples, ncol = p)
  nu_samples <- matrix(NA, nrow = num_samples, ncol = p)
  lambda_samples <- matrix(NA, nrow = num_samples, ncol = p)
  sigmaSq_samples <- numeric(num_samples)
  phi_samples <- numeric(num_samples)
  
  ### Initial values
  sigmaSq_cur <- var(y)
  phi_cur <- if (is.null(phi_val)) scale_phi else phi_val
  # beta_cur <- rnorm(p, 0, sigma = sqrt(sigmaSq_cur))
  nu_cur <- rgamma(p, 1, 1)
  lambda_cur <- rgamma(p, 1, 1)
  
  idx_store <- 0
  
  for (i in 1:num_total) {
    beta_cur <- Gibbs_beta(
      X = X, y = y, sigmaSq = sigmaSq_cur,
      phi = phi_cur, nu = nu_cur, lambda = lambda_cur,
      woodbury = woodbury, diagX = diagX
    )

    for (d in 1:p){
      chi <- max(exp(2 * log(abs(beta_cur[d])) + log(lambda_cur[d]) - log(phi_cur)),
                 10^(-306)) # 307 is smallest integer that didn't result in the code breaking here
      nu_cur[d] <- GIGrvg::rgig(1, lambda = a-0.5, chi = chi, psi = 2)
    }
    rate_lam <- 1 + beta_cur^2 / (2*phi_cur*nu_cur)
    lambda_cur <- rgamma(p, shape = b+0.5, rate  = rate_lam)
    
    sigmaSq_cur <- Gibbs_sigmaSq(n, X, y, beta_cur)
    if (!fix_phi) {
      phi_cur <- Gibbs_phi(phi_cur, beta_cur, nu_cur, lambda_cur, scale_phi)
    }
    
    # store after burn-in
    if(i > num_burnin) {
      idx_store <- idx_store + 1
      beta_samples[idx_store, ] <- beta_cur
      lambda_samples[idx_store, ] <- lambda_cur
      nu_samples[idx_store, ] <- nu_cur
      sigmaSq_samples[idx_store] <- sigmaSq_cur
      phi_samples[idx_store] <- phi_cur
    }
  }
  return(
    list(beta = beta_samples,
         lambda = lambda_samples,
         nu = nu_samples,
         sigmaSq = sigmaSq_samples,
         phi = phi_samples)
  )
}
