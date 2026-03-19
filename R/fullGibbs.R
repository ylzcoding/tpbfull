#' Fully Gibbs Sampler Function with Miller's Approximation
#'
#' @import mvtnorm
#' @import coda
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param num_output Integer, number of posterior samples to save (after burn-in)
#' @param num_burnin Integer, number of burn-in iterations to discard
#' @param woodbury Logical, use Woodbury identity in beta update
#' @param diagX Logical, assume diagonal X
#' @param hyper_prior string, "gamma" or "hcauchy"
#' @param method string, "miller" or "mh" (only applicabl. e if prior_option="gamma")
#' @param hyper_params List of hyperparameters (s_a, r_a, s_b, r_b, scale_a, scale_b)
#' @param mh_step scalar, step size for MH uniform proposal
#' @return A list containing posterior samples matrices and a diagnostics summary table.
#' @export
fullGibbs <- function(X, y, num_output = 10000, num_burnin = 10000, thin = 1,
                      woodbury = TRUE, diagX = FALSE, hyper_prior = "gamma", method = "mh", mh_step = 0.05,
                      hyper_params = list(s_a=1.5, r_a=1, s_b=1.5, r_b=1, scale_a = 1, scale_b = 1)) {

  n <- nrow(X)
  p <- ncol(X)
  total_iter <- num_output + num_burnin
  n_save <- floor(num_output / thin)

  # Initialize Parameters (Starting Values)
  sigmaSq <- var(y)
  beta    <- t(mvtnorm::rmvnorm(1, sigma = sigmaSq * diag(p)))
  a       <- 0.5
  b       <- 0.5
  phi     <- 1
  w       <- 1
  nu      <- 1
  lambda <- 1
  xi      <- rep(1, p)

  # Storage matrices (only for saved samples)
  store_beta    <- matrix(0, nrow = n_save, ncol = p)
  store_scalars <- matrix(0, nrow = n_save, ncol = 5) # sigmaSq, phi, w, a, b
  colnames(store_scalars) <- c("sigmaSq", "phi", "w", "a", "b")
  idx <- 0

  # Progress bar
  #pb <- txtProgressBar(min = 0, max = total_iter, style = 3)

  cat("Starting Gibbs Sampler...\n")

  accept_count_a <- 0
  accept_count_b <- 0

  for (iter in 1:total_iter) {
    beta <- Gibbs_beta(X = X, y = y, a = a, b = b,
                       phi = phi, sigmaSq = sigmaSq, nu = nu, lambda = lambda,
                       woodbury = woodbury, diagX = diagX)

    sigmaSq <- Gibbs_sigmaSq(n, X, y, beta)
    lambda <- Gibbs_lambda(b, p, phi, beta, nu)
    w <- Gibbs_w(phi)
    phi <- Gibbs_phi(p, nu, lambda, w, beta)
    #xi <- Gibbs_xi(a, p, nu)
    nu <- Gibbs_nu(phi, beta, lambda, a)
    a_new <- Gibbs_a(nu = nu, a = a, prior = hyper_prior,
                 s_a = hyper_params$s_a, r_a = hyper_params$r_a, scale_a = hyper_params$scale_a,
                 method = method, mh_step = mh_step)
    if (a_new != a) {
      accept_count_a <- accept_count_a + 1
    }
    a <- a_new

    b_new <- Gibbs_b(lambda = lambda, b = b, prior = hyper_prior,
                 s_b = hyper_params$s_b, r_b = hyper_params$r_b, scale_b = hyper_params$scale_b,
                 method = method, mh_step = mh_step)
    if (b_new != b) {
      accept_count_b <- accept_count_b + 1
    }
    b <- b_new

    if (iter > num_burnin) {
      if ((iter - num_burnin) %% thin == 0) {
        idx <- idx + 1
        if (idx <= n_save) {
          store_beta[idx, ] <- beta
          store_scalars[idx, ] <- c(sigmaSq, phi, w, a, b)
        }
      }
    }
    # Update progress
    # if (iter %% 100 == 0) setTxtProgressBar(pb, iter)
  }

  accept_a <- accept_count_a / total_iter
  accept_b <- accept_count_b / total_iter
  # close(pb)

  ### Diagnostics & Summary

  scalar_names <- colnames(store_scalars)
  result <- list(
    samples = list(
      beta = store_beta,
      scalars = store_scalars
    )
  )

  if (method == "mh") {
    result$acceptance_rates <- list(
      a = accept_a,
      b = accept_b
    )
  }

  return(result)
}
