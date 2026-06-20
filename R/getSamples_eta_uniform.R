#' Fully Bayesian TPB Sampler with Local eta_j ~ Uniform(0, 1)
#'
#' Baseline sampler inspired by the Yang-style TPB parameterization
#' a_j = eta_j and b_j = 1 - eta_j. This keeps the FPB extra shape parameter
#' fixed at sigma_j = 1, so the sampler remains a three-parameter TPB model
#' with local shape learning rather than the full four-parameter beta prior.
#'
#' @param X n by p design matrix.
#' @param y Response vector of length n.
#' @param num_burnin Number of burn-in iterations.
#' @param num_samples Number of posterior samples to save.
#' @param woodbury Logical; use Woodbury identity in beta updates.
#' @param diagX Logical; treat X as diagonal.
#' @param fix_phi Logical; if TRUE, phi is fixed at phi_val.
#' @param phi_val Fixed or initial value for phi.
#' @param scale_phi Half-Cauchy scale for sqrt(phi); defaults to 1.
#' @param eta_init Initial eta value. Either scalar or vector of length p.
#' @param mh_step_eta Random-walk proposal sd on logit(eta).
#' @param seed Optional random seed.
#' @return A list containing posterior samples for beta, eta, a, b, lambda,
#'   nu, sigmaSq, phi, and eta acceptance rates.
#' @export
tpb_fullyBayes_eta_uniform <- function(X, y,
                                       num_burnin,
                                       num_samples,
                                       woodbury = TRUE,
                                       diagX = FALSE,
                                       fix_phi = FALSE,
                                       phi_val = NULL,
                                       scale_phi = 1,
                                       eta_init = 0.5,
                                       mh_step_eta = 2.5,
                                       seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n <- nrow(X)
  p <- ncol(X)
  num_total <- num_burnin + num_samples
  if (num_total < 1) {
    stop("num_burnin + num_samples must be positive.")
  }
  if (!is.finite(scale_phi) || scale_phi <= 0) {
    scale_phi <- 1
  }
  if (!is.finite(mh_step_eta) || mh_step_eta <= 0) {
    stop("mh_step_eta must be positive.")
  }

  eps <- 1e-6
  eta_cur <- if (length(eta_init) == 1) {
    rep(eta_init, p)
  } else {
    eta_init
  }
  if (length(eta_cur) != p) {
    stop("eta_init must be a scalar or a vector of length ncol(X).")
  }
  eta_cur <- pmin(pmax(as.numeric(eta_cur), eps), 1 - eps)

  beta_samples <- matrix(NA_real_, nrow = num_samples, ncol = p)
  eta_samples <- matrix(NA_real_, nrow = num_samples, ncol = p)
  nu_samples <- matrix(NA_real_, nrow = num_samples, ncol = p)
  lambda_samples <- matrix(NA_real_, nrow = num_samples, ncol = p)
  sigmaSq_samples <- numeric(num_samples)
  phi_samples <- numeric(num_samples)

  sigmaSq_cur <- max(stats::var(as.vector(y)), .Machine$double.xmin)
  phi_cur <- if (is.null(phi_val)) scale_phi else phi_val
  phi_cur <- max(phi_cur, .Machine$double.xmin)
  beta_cur <- rep(0, p)
  nu_cur <- stats::rgamma(p, shape = 1, rate = 1)
  lambda_cur <- stats::rgamma(p, shape = 1, rate = 1)

  eta_accept <- numeric(p)
  idx_store <- 0

  log_eta_kernel <- function(eta, nu, lambda) {
    if (!is.finite(eta) || eta <= 0 || eta >= 1) {
      return(-Inf)
    }
    stats::dgamma(nu, shape = eta, rate = 1, log = TRUE) +
      stats::dgamma(lambda, shape = 1 - eta, rate = 1, log = TRUE)
  }

  log_eta_z_kernel <- function(z, nu, lambda) {
    eta <- stats::plogis(z)
    log_eta_kernel(eta, nu, lambda) + log(eta) + log1p(-eta)
  }

  fitted_values <- function(beta) {
    if (diagX) {
      diag(X) * beta
    } else {
      as.vector(X %*% beta)
    }
  }

  for (iter in seq_len(num_total)) {
    beta_cur <- Gibbs_beta(
      X = X, y = y,
      phi = phi_cur,
      sigmaSq = sigmaSq_cur,
      nu = nu_cur,
      lambda = lambda_cur,
      woodbury = woodbury,
      diagX = diagX
    )

    for (j in seq_len(p)) {
      z_cur <- stats::qlogis(eta_cur[j])
      z_prop <- z_cur + stats::rnorm(1, sd = mh_step_eta)
      log_acc <- log_eta_z_kernel(z_prop, nu_cur[j], lambda_cur[j]) -
        log_eta_z_kernel(z_cur, nu_cur[j], lambda_cur[j])
      if (is.finite(log_acc) && log(stats::runif(1)) < log_acc) {
        eta_cur[j] <- stats::plogis(z_prop)
        eta_accept[j] <- eta_accept[j] + 1
      }
    }
    eta_cur <- pmin(pmax(eta_cur, eps), 1 - eps)

    chi <- pmax(beta_cur^2 * lambda_cur / phi_cur, .Machine$double.xmin)
    nu_cur <- vapply(seq_len(p), function(j) {
      GIGrvg::rgig(1, lambda = eta_cur[j] - 0.5, chi = chi[j], psi = 2)
    }, numeric(1))
    nu_cur <- pmax(nu_cur, .Machine$double.xmin)

    lambda_rate <- 1 + beta_cur^2 / (2 * phi_cur * nu_cur)
    lambda_cur <- stats::rgamma(
      p,
      shape = 1 - eta_cur + 0.5,
      rate = pmax(lambda_rate, .Machine$double.xmin)
    )
    lambda_cur <- pmax(lambda_cur, .Machine$double.xmin)

    resid <- as.vector(y) - fitted_values(beta_cur)
    sigmaSq_cur <- 1 / stats::rgamma(
      1,
      shape = n / 2,
      rate = max(sum(resid^2) / 2, .Machine$double.xmin)
    )
    sigmaSq_cur <- max(sigmaSq_cur, .Machine$double.xmin)

    if (!fix_phi) {
      w <- 1 / stats::rgamma(
        1,
        shape = 1,
        rate = 1 / scale_phi + 1 / phi_cur
      )
      phi_rate <- 1 / w + sum(lambda_cur * beta_cur^2 / (2 * nu_cur))
      phi_cur <- 1 / stats::rgamma(
        1,
        shape = (p + 1) / 2,
        rate = max(phi_rate, .Machine$double.xmin)
      )
      phi_cur <- max(phi_cur, .Machine$double.xmin)
    }

    if (iter > num_burnin) {
      idx_store <- idx_store + 1
      beta_samples[idx_store, ] <- beta_cur
      eta_samples[idx_store, ] <- eta_cur
      nu_samples[idx_store, ] <- nu_cur
      lambda_samples[idx_store, ] <- lambda_cur
      sigmaSq_samples[idx_store] <- sigmaSq_cur
      phi_samples[idx_store] <- phi_cur
    }
  }

  list(
    beta = beta_samples,
    eta = eta_samples,
    a = eta_samples,
    b = 1 - eta_samples,
    lambda = lambda_samples,
    nu = nu_samples,
    sigmaSq = sigmaSq_samples,
    phi = phi_samples,
    eta_acceptance = eta_accept / num_total,
    sigma_i = rep(1, p),
    mh_step_eta = mh_step_eta
  )
}
