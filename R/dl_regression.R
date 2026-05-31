#' Fit Bayesian linear regression with the Dirichlet-Laplace prior
#'
#' Gibbs sampler for the Dirichlet-Laplace shrinkage prior of Bhattacharya
#' et al. (2015).  The regression model is y = X beta + e with no intercept;
#' callers should center y and standardize X when that is desired.
#'
#' @param X n by p design matrix.
#' @param y Response vector of length n.
#' @param num_burnin Number of warmup iterations.
#' @param num_samples Number of post-warmup iterations.
#' @param thin Save every thin-th post-warmup draw.
#' @param hyper Dirichlet concentration parameter.  If NULL, use 1 / max(n, p).
#' @param sigma_a,sigma_b Shape and scale of the inverse-gamma prior on sigma^2.
#'   The default IG(1.5, 0.5) matches the variance prior used in the other
#'   baseline samplers in these simulations.
#' @param woodbury Use the fast p > n Gaussian update.
#' @param seed Optional RNG seed.
#' @return A list containing posterior beta samples and scalar traces.
#' @export
dl_regression <- function(X, y,
                          num_burnin = 5000,
                          num_samples = 5000,
                          thin = 1,
                          hyper = NULL,
                          sigma_a = 1.5,
                          sigma_b = 0.5,
                          woodbury = TRUE,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  p <- ncol(X)
  if (length(y) != n) stop("Length of y must match nrow(X).")
  if (is.null(hyper)) hyper <- 1 / max(n, p)
  if (hyper <= 0) stop("hyper must be positive.")

  total_iter <- num_burnin + num_samples
  save_iter <- seq.int(num_burnin + thin, total_iter, by = thin)
  n_save <- length(save_iter)

  beta_samples <- matrix(NA_real_, nrow = n_save, ncol = p)
  sigma2_samples <- numeric(n_save)
  tau_samples <- numeric(n_save)

  beta <- rep(0, p)
  sigma2 <- max(stats::var(y), .Machine$double.eps)
  local_psi <- rep(1, p)
  phi <- rep(1 / p, p)
  tau <- 1

  xtx <- crossprod(X)
  xty <- crossprod(X, y)
  save_id <- 0L

  for (iter in seq_len(total_iter)) {
    local_var <- pmax(local_psi * phi^2 * tau^2, .Machine$double.xmin)

    beta <- dl_draw_beta(
      X = X, y = y, xtx = xtx, xty = xty,
      sigma2 = sigma2, local_var = local_var,
      woodbury = woodbury
    )

    resid <- y - drop(X %*% beta)
    quad_prior <- sum(beta^2 / local_var)
    sigma2 <- 1 / stats::rgamma(
      1,
      shape = sigma_a + 0.5 * (n + p),
      rate = sigma_b + 0.5 * (sum(resid^2) + quad_prior)
    )
    sigma <- sqrt(sigma2)

    abs_beta <- pmax(abs(beta), .Machine$double.xmin)

    inv_psi_mean <- pmax(sigma * phi * tau / abs_beta, .Machine$double.xmin)
    inv_psi <- statmod::rinvgauss(p, mean = inv_psi_mean, shape = 1)
    local_psi <- 1 / pmax(inv_psi, .Machine$double.xmin)

    # R2D2::dl samples tau | phi, beta with chi = 2 sum_j |beta_j| /
    # (sigma phi_j), followed by phi | beta with chi_j = 2 |beta_j| / sigma.
    # The phi conditional does not include tau after normalizing the latent
    # T_j variables to the simplex.
    tau_chi <- pmax(2 * sum(abs_beta / (sigma * pmax(phi, .Machine$double.xmin))),
                    .Machine$double.xmin)
    tau <- GIGrvg::rgig(
      1,
      lambda = p * (hyper - 1),
      chi = tau_chi,
      psi = 1
    )
    tau <- max(tau, .Machine$double.xmin)

    phi_gig_chi <- pmax(2 * abs_beta / sigma, .Machine$double.xmin)
    phi_raw <- vapply(
      phi_gig_chi,
      function(chi) GIGrvg::rgig(1, lambda = hyper - 1, chi = chi, psi = 1),
      numeric(1)
    )
    phi <- phi_raw / sum(phi_raw)

    if (iter %in% save_iter) {
      save_id <- save_id + 1L
      beta_samples[save_id, ] <- beta
      sigma2_samples[save_id] <- sigma2
      tau_samples[save_id] <- tau
    }
  }

  colnames(beta_samples) <- colnames(X)
  list(
    beta = beta_samples,
    sigma2_samples = sigma2_samples,
    tau_samples = tau_samples,
    hyper = hyper
  )
}

dl_draw_beta <- function(X, y, xtx, xty, sigma2, local_var, woodbury = TRUE) {
  n <- nrow(X)
  p <- ncol(X)

  if (woodbury && p > n) {
    u <- stats::rnorm(p, 0, sqrt(sigma2 * local_var))
    delta <- stats::rnorm(n)
    XD <- sweep(X, 2, local_var, FUN = "*")
    M <- XD %*% t(X) + diag(n)
    rhs <- y - drop(X %*% u) - sqrt(sigma2) * delta
    w <- solve(M, rhs)
    return(u + local_var * drop(crossprod(X, w)))
  }

  precision <- xtx + diag(1 / local_var, p)
  chol_precision <- chol(precision)
  mean_beta <- backsolve(chol_precision, forwardsolve(t(chol_precision), xty))
  z <- stats::rnorm(p)
  drop(mean_beta + sqrt(sigma2) * backsolve(chol_precision, z))
}
