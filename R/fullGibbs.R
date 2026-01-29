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
#' @param hyper_params List of hyperparameters (s_a, r_a, s_b, r_b, sigma_a, sigma_b)
#' @return A list containing posterior samples matrices and a diagnostics summary table.
#' @export
fullGibbs <- function(X, y, num_output = 10000, num_burnin = 10000,
                      woodbury = FALSE, diagX = FALSE,
                      hyper_params = list(s_a=1, r_a=1, s_b=1, r_b=1)) {

  n <- nrow(X)
  p <- ncol(X)
  total_iter <- num_output + num_burnin

  # Initialize Parameters (Starting Values)
  sigmaSq <- var(y)
  beta    <- t(mvtnorm::rmvnorm(1, sigma = sigmaSq * diag(p)))
  a       <- 0.5
  b       <- 0.5
  phi     <- 1
  w       <- 1
  nu      <- rgamma(p, a, 1)
  lambda <- rgamma(p, b, 1)
  xi      <- rep(1, p)

  # Storage matrices (only for saved samples)
  store_beta    <- matrix(0, nrow = num_output, ncol = p)
  store_scalars <- matrix(0, nrow = num_output, ncol = 5) # sigmaSq, phi, w, a, b
  colnames(store_scalars) <- c("sigmaSq", "phi", "w", "a", "b")

  # Progress bar
  #pb <- txtProgressBar(min = 0, max = total_iter, style = 3)

  cat("Starting Gibbs Sampler...\n")
  cat(sprintf("  N: %d, P: %d\n", n, p))
  cat(sprintf("  Burn-in: %d, Saved: %d\n", num_burnin, num_output))

  for (iter in 1:total_iter) {
    beta <- Gibbs_beta(X = X, y = y, a = a, b = b,
                       phi = phi, sigmaSq = sigmaSq, nu = nu, lambda = lambda,
                       woodbury = woodbury, diagX = diagX)

    sigmaSq <- Gibbs_sigmaSq(n, X, y, beta)
    lambda <- Gibbs_lambda(b, p, phi, beta, nu)
    if (any(is.na(lambda)) || any(lambda <= 0)) {
      warning(sprintf("Iteration %d: 'lambda' contains NA or non-positive values. Replacing <=0 with 1e-10 for stability.", iter))
      lambda[is.na(lambda)] <- 1e-10
      lambda[lambda <= 0] <- 1e-10
    }

    w <- Gibbs_w(phi)
    phi <- Gibbs_phi(p, nu, lambda, w, beta)
    #xi <- Gibbs_xi(a, p, nu)
    nu <- Gibbs_nu(phi, beta, lambda, a)
    if (any(is.na(nu)) || any(nu <= 0)) {
      # Try to fix numerically unstable values (optional) or stop to debug
      warning(sprintf("Iteration %d: 'nu' contains NA or non-positive values. Replacing <=0 with 1e-10 for stability.", iter))
      nu[is.na(nu)] <- 1e-10
      nu[nu <= 0] <- 1e-10
      # If you prefer to stop debugging immediately, uncomment the line below:
      # stop("Critical Error: 'nu' became invalid.")
    }

    a <- Gibbs_a(nu = nu, a = a, s_a = hyper_params$s_a, r_a = hyper_params$r_a)
    b <- Gibbs_b(lambda = lambda, b = b, s_b = hyper_params$s_b, r_b = hyper_params$r_b)

    if (iter > num_burnin) {
      idx <- iter - num_burnin
      store_beta[idx, ] <- beta
      store_scalars[idx, ] <- c(sigmaSq, phi, w, a, b)
    }

    # Update progress
    # if (iter %% 100 == 0) setTxtProgressBar(pb, iter)
  }

  # close(pb)

  ### Diagnostics & Summary

  # Combine scalars for easier handling
  scalar_names <- colnames(store_scalars)

  return(list(
    samples = list(
      beta = store_beta,
      scalars = store_scalars
    )
  ))
}
