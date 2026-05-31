eb_positive_finite <- function(x, lower) {
  if (!is.finite(x) || x <= lower) {
    return(lower)
  }
  x
}

eb_isConverged_hyper <- function(..., delta1, delta2, window_size = 3) {
  vecs <- list(...)
  if (any(vapply(vecs, length, integer(1)) < window_size + 1)) {
    return(FALSE)
  }
  rel_err <- function(vec) {
    w <- tail(vec, window_size + 1)
    rel <- abs(diff(w)) / (abs(head(w, -1)) + delta1)
    if (all(is.na(rel))) Inf else max(rel, na.rm = TRUE)
  }
  max(vapply(vecs, rel_err, numeric(1)), na.rm = TRUE) < delta2
}

eb_olasso_screen <- function(X, y, epsilon = 1e-8) {
  ol_fit <- natural::olasso(X, y, intercept = FALSE)
  list(
    beta = ol_fit$beta_1,
    sigmaSq = max(epsilon, ol_fit$sig_obj_1),
    source = "olasso"
  )
}

#' Default Candidate Set for TPB Model Competition
#'
#' @return A named list of candidate (a,b) pairs.
#' @export
tpb_default_candidates <- function() {
  list(
    hs = list(a = 0.5, b = 0.5),
    lasso = list(a = 1.0, b = 5.0),
    normal_gamma = list(a = 0.5, b = 5.0),
    ridge = list(a = 5.0, b = 5.0)
  )
}

eb_ridge_screen <- function(X, y, epsilon = 1e-8, ridge_lambda = NULL,
                            diagX = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  y <- as.vector(y)
  if (is.null(ridge_lambda)) {
    ridge_lambda <- sqrt(n)
  }
  ridge_lambda <- eb_positive_finite(ridge_lambda, epsilon)

  if (diagX) {
    x_diag <- diag(X)
    beta <- x_diag * y / (x_diag^2 + ridge_lambda)
    fitted <- x_diag * beta
    df <- sum(x_diag^2 / (x_diag^2 + ridge_lambda))
  } else if (p <= n) {
    beta <- solve(crossprod(X) + ridge_lambda * diag(p), crossprod(X, y))
    fitted <- as.vector(X %*% beta)
    singular_values <- svd(X, nu = 0, nv = 0)$d^2
    df <- sum(singular_values / (singular_values + ridge_lambda))
  } else {
    alpha <- solve(tcrossprod(X) + ridge_lambda * diag(n), y)
    beta <- as.vector(crossprod(X, alpha))
    fitted <- as.vector(X %*% beta)
    singular_values <- svd(X, nu = 0, nv = 0)$d^2
    df <- sum(singular_values / (singular_values + ridge_lambda))
  }

  rss <- sum((y - fitted)^2)
  sigmaSq <- rss / max(n - df, 1)
  list(
    beta = as.vector(beta),
    sigmaSq = max(epsilon, sigmaSq),
    source = "ridge",
    ridge_lambda = ridge_lambda
  )
}

eb_calculate_marginal_loglik_beta <- function(beta_vec, model_params) {
  a <- model_params$a
  b <- model_params$b
  # EB stores omega in beta | nu, lambda, omega.
  # The marginal TPB beta density uses phi = omega * b / a.
  phi <- model_params$omega * b / a
  if (!is.finite(phi) || phi <= 0) {
    return(-Inf)
  }
  dens <- eb_d_tpb(beta_vec = beta_vec, a = a, b = b, phi = phi)
  if (any(is.na(dens)) || any(dens <= 0)) {
    return(-Inf)
  }
  loglik_individual <- log(pmin(pmax(dens, .Machine$double.xmin),
                                .Machine$double.xmax))
  sum(loglik_individual)
}

eb_calculate_gaussian_loglik_beta <- function(beta_vec, X, y, sigmaSq,
                                             diagX = FALSE) {
  sigmaSq <- eb_positive_finite(sigmaSq, 1e-16)
  fitted <- if (diagX) {
    diag(X) * beta_vec
  } else {
    X %*% beta_vec
  }
  resid <- as.vector(y) - as.vector(fitted)
  -0.5 * length(resid) * log(2 * pi * sigmaSq) -
    sum(resid^2) / (2 * sigmaSq)
}

eb_calculate_posterior_kernel_score <- function(beta_vec, model_params, X, y,
                                                diagX = FALSE) {
  eb_calculate_gaussian_loglik_beta(
    beta_vec = beta_vec,
    X = X,
    y = y,
    sigmaSq = model_params$sigmaSq,
    diagX = diagX
  ) +
    eb_calculate_marginal_loglik_beta(
      beta_vec = beta_vec,
      model_params = model_params
    )
}

eb_logsumexp <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(-Inf)
  }
  m <- max(x)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(x - m)))
}

eb_importance_ess <- function(logw) {
  logw <- logw[!is.na(logw)]
  if (length(logw) == 0) {
    return(0)
  }
  log_sum_w <- eb_logsumexp(logw)
  log_sum_w2 <- eb_logsumexp(2 * logw)
  ess <- exp(2 * log_sum_w - log_sum_w2)
  if (is.finite(ess)) ess else 0
}

eb_complete_logpost <- function(a, b, omega, sigmaSq, X, y,
                                beta, nu, lambda, xi = NULL,
                                diagX = FALSE) {
  a <- eb_positive_finite(a, 1e-6)
  b <- eb_positive_finite(b, 1e-6)
  omega <- eb_positive_finite(omega, 1e-16)
  sigmaSq <- eb_positive_finite(sigmaSq, 1e-16)
  nu <- pmax(nu, .Machine$double.xmin)
  lambda <- pmax(lambda, .Machine$double.xmin)

  local_sd <- sqrt(pmax(omega * nu / lambda, .Machine$double.xmin))
  log_xi <- if (is.null(xi)) {
    0
  } else {
    xi <- pmax(xi, .Machine$double.xmin)
    sum(stats::dgamma(xi, shape = a, rate = 1 / nu, log = TRUE))
  }

  eb_calculate_gaussian_loglik_beta(beta, X, y, sigmaSq, diagX) +
    sum(stats::dnorm(beta, mean = 0, sd = local_sd, log = TRUE)) +
    sum(stats::dgamma(nu, shape = a, rate = a, log = TRUE)) +
    sum(stats::dgamma(lambda, shape = b, rate = b, log = TRUE)) +
    log_xi
}

eb_em_importance_weights <- function(samples, target_params, source_params,
                                     X, y, diagX = FALSE) {
  n_samples <- nrow(samples$beta)
  log_target <- vapply(seq_len(n_samples), function(i) {
    eb_complete_logpost(
      a = target_params$a, b = target_params$b,
      omega = target_params$omega, sigmaSq = target_params$sigmaSq,
      X = X, y = y,
      beta = samples$beta[i, ],
      nu = samples$nu[i, ],
      lambda = samples$lambda[i, ],
      xi = if (is.null(samples$xi)) NULL else samples$xi[i, ],
      diagX = diagX
    )
  }, numeric(1))
  log_source <- vapply(seq_len(n_samples), function(i) {
    eb_complete_logpost(
      a = source_params$a, b = source_params$b,
      omega = source_params$omega, sigmaSq = source_params$sigmaSq,
      X = X, y = y,
      beta = samples$beta[i, ],
      nu = samples$nu[i, ],
      lambda = samples$lambda[i, ],
      xi = if (is.null(samples$xi)) NULL else samples$xi[i, ],
      diagX = diagX
    )
  }, numeric(1))

  logw <- log_target - log_source
  ess <- eb_importance_ess(logw)
  if (all(is.na(logw)) || all(!is.finite(logw))) {
    weights <- rep(1, n_samples)
    ess <- 0
  } else {
    max_logw <- max(logw, na.rm = TRUE)
    weights <- exp(logw - max_logw)
    weights[!is.finite(weights)] <- 0
    if (!is.finite(sum(weights)) || sum(weights) <= 0) {
      weights <- rep(1, n_samples)
      ess <- 0
    }
  }

  list(weights = weights, ess = ess, log_weights = logw)
}

eb_weighted_m_step_sigmaSq <- function(beta_matrix, weights, X, y,
                                       diagX = FALSE) {
  weights <- as.numeric(weights)
  w_sum <- sum(weights)
  if (!is.finite(w_sum) || w_sum <= 0) {
    return(eb_m_step_sigmaSq(beta_matrix, X, y, diagX))
  }

  y <- as.vector(y)
  n <- nrow(X)
  num_samples <- nrow(beta_matrix)
  ssr <- if (diagX) {
    x_diag <- diag(X)
    vapply(seq_len(num_samples), function(i) {
      sum((y - x_diag * beta_matrix[i, ])^2)
    }, numeric(1))
  } else {
    vapply(seq_len(num_samples), function(i) {
      sum((y - X %*% beta_matrix[i, ])^2)
    }, numeric(1))
  }

  eb_positive_finite(sum(ssr * weights) / (w_sum * n), 1e-16)
}

eb_weighted_m_step_omega <- function(beta_matrix, lambda_matrix, nu_matrix,
                                     weights) {
  weights <- as.numeric(weights)
  vals <- beta_matrix^2 * lambda_matrix / pmax(nu_matrix, 1e-16)
  row_vals <- rowMeans(vals, na.rm = TRUE)
  ok <- is.finite(row_vals) & is.finite(weights) & weights >= 0
  denom <- sum(weights[ok])
  if (!is.finite(denom) || denom <= 0) {
    return(eb_m_step_omega(beta_matrix, lambda_matrix, nu_matrix))
  }
  eb_positive_finite(sum(row_vals[ok] * weights[ok]) / denom, 1e-16)
}

eb_softmax <- function(log_weights) {
  if (all(!is.finite(log_weights))) {
    return(rep(1 / length(log_weights), length(log_weights)))
  }
  m <- max(log_weights, na.rm = TRUE)
  w <- exp(log_weights - m)
  w[!is.finite(w)] <- 0
  sw <- sum(w)
  if (is.finite(sw) && sw > 0) {
    w / sw
  } else {
    rep(1 / length(log_weights), length(log_weights))
  }
}

eb_precompute_selection_scores <- function(presampled_betas,
                                           pre_optimized_params,
                                           X, y,
                                           diagX = FALSE) {
  model_names <- names(pre_optimized_params)
  scores <- vector("list", length(model_names))
  names(scores) <- model_names

  for (source_name in model_names) {
    beta_mat <- presampled_betas[[source_name]]
    n_s <- nrow(beta_mat)
    score_mat <- matrix(NA_real_, nrow = n_s, ncol = length(model_names),
                        dimnames = list(NULL, model_names))
    for (target_name in model_names) {
      score_mat[, target_name] <- vapply(seq_len(n_s), function(i) {
        eb_calculate_posterior_kernel_score(
          beta_vec = beta_mat[i, ],
          model_params = pre_optimized_params[[target_name]],
          X = X,
          y = y,
          diagX = diagX
        )
      }, numeric(1))
    }
    scores[[source_name]] <- score_mat
  }

  scores
}

eb_reverse_logistic_selection <- function(presampled_betas,
                                          pre_optimized_params,
                                          X, y,
                                          diagX = FALSE,
                                          kernel_scores = NULL) {
  model_names <- names(pre_optimized_params)
  M <- length(model_names)
  if (M == 1) {
    return(list(
      model_probabilities = stats::setNames(1, model_names),
      log_evidence = stats::setNames(0, model_names),
      convergence = 0,
      objective = 0,
      method = "importance_reverse_logistic"
    ))
  }

  if (is.null(kernel_scores)) {
    kernel_scores <- eb_precompute_selection_scores(
      presampled_betas = presampled_betas,
      pre_optimized_params = pre_optimized_params,
      X = X,
      y = y,
      diagX = diagX
    )
  }

  n_by_source <- vapply(presampled_betas, nrow, integer(1))
  log_n <- log(pmax(n_by_source, 1))
  ref <- 1

  objective_gradient <- function(par) {
    theta <- numeric(M)
    theta[-ref] <- par
    value <- 0
    grad <- numeric(M)

    for (source_idx in seq_len(M)) {
      score_mat <- kernel_scores[[model_names[source_idx]]]
      adjusted <- sweep(score_mat, 2, theta, FUN = "-")
      adjusted <- sweep(adjusted, 2, log_n, FUN = "+")

      for (i in seq_len(nrow(adjusted))) {
        row <- adjusted[i, ]
        log_den <- eb_logsumexp(row)
        if (!is.finite(log_den)) {
          value <- value + 1e8
          next
        }
        probs <- exp(row - log_den)
        probs[!is.finite(probs)] <- 0
        value <- value + theta[source_idx] + log_den
        grad[source_idx] <- grad[source_idx] + 1
        grad <- grad - probs
      }
    }

    list(value = value, gradient = grad[-ref])
  }

  opt <- tryCatch({
    stats::optim(
      par = rep(0, M - 1),
      fn = function(par) objective_gradient(par)$value,
      gr = function(par) objective_gradient(par)$gradient,
      method = "BFGS",
      control = list(maxit = 1000, reltol = 1e-8)
    )
  }, error = function(e) NULL)

  theta <- numeric(M)
  convergence <- NA_integer_
  objective <- NA_real_
  if (!is.null(opt) && all(is.finite(opt$par))) {
    theta[-ref] <- opt$par
    convergence <- opt$convergence
    objective <- opt$value
  }

  if (any(!is.finite(theta))) {
    theta <- rep(0, M)
  }
  theta <- theta - max(theta, na.rm = TRUE)
  probs <- eb_softmax(theta)
  names(theta) <- model_names
  names(probs) <- model_names

  list(
    model_probabilities = probs,
    log_evidence = theta,
    n_by_source = n_by_source,
    convergence = convergence,
    objective = objective,
    method = "importance_reverse_logistic"
  )
}

eb_d_tpb <- function(beta_vec, a, b, phi) {
  if (!is.finite(a) || !is.finite(b) || !is.finite(phi) ||
      a <= 0 || b <= 0 || phi <= 0) {
    return(rep(NA_real_, length(beta_vec)))
  }

  tricomi_U_integral <- function(U_a, U_b, z) {
    integrand <- function(t) {
      exp(-z * t) * t^(U_a - 1) * (1 + t)^(U_b - U_a - 1)
    }
    integral_result <- tryCatch({
      integrate(integrand, lower = 0, upper = Inf)$value
    }, error = function(e) NA_real_)
    integral_result / gamma(U_a)
  }

  log_const <- lgamma(0.5 + b) + lgamma(a + b) -
    lgamma(a) - lgamma(b) - 0.5 * log(2 * pi * phi)
  U_a <- 0.5 + b
  U_b <- 1.5 - a

  z_vals <- beta_vec^2 / (2 * phi)
  out <- rep(NA_real_, length(beta_vec))
  valid <- is.finite(z_vals) & z_vals >= 0
  if (!any(valid)) {
    return(out)
  }

  U_vals <- rep(NA_real_, length(beta_vec))
  U_vals[valid] <- tryCatch({
    gsl::hyperg_U(U_a, U_b, z_vals[valid])
  }, error = function(e) rep(NA_real_, sum(valid)))

  needs_fallback <- valid & (!is.finite(U_vals) | U_vals <= 0)
  if (any(needs_fallback)) {
    U_vals[needs_fallback] <- vapply(z_vals[needs_fallback], function(z_val) {
      tricomi_U_integral(U_a, U_b, z_val)
    }, numeric(1))
  }

  ok <- valid & is.finite(U_vals) & U_vals > 0
  out[ok] <- exp(pmin(log_const + log(U_vals[ok]), log(.Machine$double.xmax)))
  out
}

eb_getsamples_emp <- function(num, X, y, a, b, omega, sigmaSq,
                              beta0, nu0, lambda0, xi0, burn = 0,
                              woodbury = FALSE, diagX = FALSE) {
  total <- num + burn
  if (total < 1) {
    stop("num + burn must be at least 1.")
  }

  p <- ncol(X)
  beta_samples <- matrix(NA_real_, nrow = total, ncol = p)
  nu_samples <- matrix(NA_real_, nrow = total, ncol = p)
  lambda_samples <- matrix(NA_real_, nrow = total, ncol = p)
  xi_samples <- matrix(NA_real_, nrow = total, ncol = p)

  beta_samples[1, ] <- beta0
  nu_samples[1, ] <- pmax(nu0, 1e-16)
  lambda_samples[1, ] <- pmax(lambda0, 1e-16)
  xi_samples[1, ] <- pmax(xi0, 1e-16)

  if (total >= 2) {
    for (i in 2:total) {
      omega_current <- eb_positive_finite(omega, 1e-16)
      beta_samples[i, ] <- Gibbs_beta(
        X = X, y = y, a = a, b = b,
        phi = omega_current,
        sigmaSq = eb_positive_finite(sigmaSq, 1e-16),
        nu = nu_samples[i - 1, ],
        lambda = lambda_samples[i - 1, ],
        woodbury = woodbury, diagX = diagX
      )
      nu_samples[i, ] <- Gibbs_nu(
        p = p, phi = omega_current, beta = beta_samples[i, ],
        lambda = lambda_samples[i - 1, ], xi = xi_samples[i - 1, ],
        a = a
      )
      xi_samples[i, ] <- Gibbs_xi(
        p = p, a = a, nu = nu_samples[i, ]
      )
      lambda_samples[i, ] <- Gibbs_lambda(
        p = p, b = b, phi = omega_current,
        beta = beta_samples[i, ], nu = nu_samples[i, ]
      )
    }
  }

  keep <- (burn + 1):total
  list(
    beta = beta_samples[keep, , drop = FALSE],
    nu = nu_samples[keep, , drop = FALSE],
    lambda = lambda_samples[keep, , drop = FALSE],
    xi = xi_samples[keep, , drop = FALSE]
  )
}

eb_m_step_sigmaSq <- function(beta_matrix, X, y, diagX = FALSE) {
  y <- as.vector(y)
  n <- nrow(X)
  num_samples <- nrow(beta_matrix)

  ssr <- if (diagX) {
    x_diag <- diag(X)
    vapply(seq_len(num_samples), function(i) {
      sum((y - x_diag * beta_matrix[i, ])^2)
    }, numeric(1))
  } else {
    vapply(seq_len(num_samples), function(i) {
      sum((y - X %*% beta_matrix[i, ])^2)
    }, numeric(1))
  }

  eb_positive_finite(mean(ssr) / n, 1e-16)
}

eb_m_step_omega <- function(beta_matrix, lambda_matrix, nu_matrix) {
  vals <- beta_matrix^2 * lambda_matrix / pmax(nu_matrix, 1e-16)
  finite_vals <- vals[is.finite(vals)]
  if (length(finite_vals) == 0) {
    return(1e-6)
  }
  eb_positive_finite(mean(finite_vals), 1e-16)
}

eb_run_em_engine <- function(X, y,
                             a_init, b_init, omega_init, sigmaSq_init,
                             max_iter,
                             iter_burnin, iter_samples,
                             delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                             window_size = 5,
                             woodbury = TRUE,
                             IS_period = 1,
                             min_em_ess_fraction = 0.1,
                             diagX = FALSE,
                             verbose = TRUE) {
  p <- ncol(X)
  IS_period <- as.integer(IS_period)
  if (!is.finite(IS_period) || IS_period < 1) {
    stop("IS_period must be a positive integer.")
  }
  if (!is.finite(min_em_ess_fraction) || min_em_ess_fraction < 0) {
    stop("min_em_ess_fraction must be a non-negative number.")
  }

  a_vec <- c(eb_positive_finite(a_init, 1e-6))
  b_vec <- c(eb_positive_finite(b_init, 1e-6))
  omega_vec <- c(eb_positive_finite(omega_init, 1e-16))
  sigmaSq_vec <- c(eb_positive_finite(sigmaSq_init, 1e-16))
  beta_diff_vec <- numeric(0)
  em_is_ess <- numeric(0)
  em_is_used <- logical(0)
  em_is_resampled <- logical(0)
  em_is_fallback <- logical(0)

  beta0 <- rnorm(p, mean = 0, sd = sqrt(sigmaSq_vec[1]))
  nu0 <- rgamma(p, 1, 1)
  lambda0 <- rgamma(p, 1, 1)
  xi0 <- Gibbs_xi(p = p, a = a_vec[1], nu = nu0)
  last_beta_hat <- NULL
  source_samples <- NULL
  source_params <- NULL

  for (k in seq_len(max_iter)) {
    current_params <- list(
      a = a_vec[k], b = b_vec[k],
      omega = omega_vec[k], sigmaSq = sigmaSq_vec[k]
    )
    should_resample <- is.null(source_samples) || ((k - 1) %% IS_period == 0)
    weights <- NULL
    used_importance <- FALSE
    fallback_to_resample <- FALSE

    if (!should_resample) {
      iw <- eb_em_importance_weights(
        samples = source_samples,
        target_params = current_params,
        source_params = source_params,
        X = X,
        y = y,
        diagX = diagX
      )
      min_ess <- max(3, min_em_ess_fraction * nrow(source_samples$beta))
      if (is.finite(iw$ess) && iw$ess >= min_ess) {
        samples <- source_samples
        weights <- iw$weights
        used_importance <- TRUE
      } else {
        should_resample <- TRUE
        fallback_to_resample <- TRUE
      }
      em_is_ess <- c(em_is_ess, iw$ess)
    }

    if (should_resample) {
      samples <- eb_getsamples_emp(
        num = iter_samples, X = X, y = y,
        a = a_vec[k], b = b_vec[k],
        omega = omega_vec[k], sigmaSq = sigmaSq_vec[k],
        beta0 = beta0, nu0 = nu0, lambda0 = lambda0, xi0 = xi0,
        burn = iter_burnin,
        woodbury = woodbury, diagX = diagX
      )
      source_samples <- samples
      source_params <- current_params
      beta0 <- samples$beta[nrow(samples$beta), ]
      nu0 <- samples$nu[nrow(samples$nu), ]
      lambda0 <- samples$lambda[nrow(samples$lambda), ]
      xi0 <- samples$xi[nrow(samples$xi), ]
      if (!fallback_to_resample) {
        em_is_ess <- c(em_is_ess, nrow(samples$beta))
      }
    }

    if (is.null(weights)) {
      sigmaSq_new <- eb_m_step_sigmaSq(samples$beta, X, y, diagX)
      omega_new <- eb_m_step_omega(samples$beta, samples$lambda, samples$nu)
      beta_hat <- colMeans(samples$beta)
    } else {
      sigmaSq_new <- eb_weighted_m_step_sigmaSq(samples$beta, weights, X, y, diagX)
      omega_new <- eb_weighted_m_step_omega(samples$beta, samples$lambda, samples$nu, weights)
      beta_hat <- colSums(sweep(samples$beta, 1, weights, FUN = "*")) / sum(weights)
    }

    em_is_used <- c(em_is_used, used_importance)
    em_is_resampled <- c(em_is_resampled, should_resample)
    em_is_fallback <- c(em_is_fallback, fallback_to_resample)
    beta_diff_vec <- c(
      beta_diff_vec,
      if (is.null(last_beta_hat)) NA_real_ else mean((beta_hat - last_beta_hat)^2)
    )
    last_beta_hat <- beta_hat

    a_vec <- c(a_vec, a_vec[k])
    b_vec <- c(b_vec, b_vec[k])
    sigmaSq_vec <- c(sigmaSq_vec, sigmaSq_new)
    omega_vec <- c(omega_vec, omega_new)

    hyper_converged <- eb_isConverged_hyper(
      a_vec, b_vec, sigmaSq_vec, omega_vec,
      delta1 = delta1, delta2 = delta2, window_size = window_size
    )
    beta_converged <- length(beta_diff_vec) >= window_size &&
      max(tail(beta_diff_vec, window_size), na.rm = TRUE) < delta3

    if (hyper_converged && beta_converged) {
      if (isTRUE(verbose)) {
        cat("Engine converged after", k, "iterations.\n")
      }
      break
    }

    if (isTRUE(verbose)) {
      cat(k, "-th iteration completed.\n")
      cat("a=", a_vec[k + 1], "\n")
      cat("b=", b_vec[k + 1], "\n")
      cat("omega=", omega_vec[k + 1], "\n")
      cat("sigmaSq=", sigmaSq_vec[k + 1], "\n")
    }
  }

  l <- length(a_vec)
  list(
    params = list(a = a_vec[l], b = b_vec[l],
                  omega = omega_vec[l], sigmaSq = sigmaSq_vec[l]),
    final_state = list(beta0 = beta0, nu0 = nu0, lambda0 = lambda0, xi0 = xi0),
    trajectories = list(a_traj = a_vec, b_traj = b_vec,
                        sigmaSq_traj = sigmaSq_vec,
                        omega_traj = omega_vec,
                        beta_diff_traj = beta_diff_vec,
                        em_is = list(
                          ess = em_is_ess,
                          used_importance = em_is_used,
                          resampled = em_is_resampled,
                          fallback_resample = em_is_fallback,
                          IS_period = IS_period,
                          min_ess_fraction = min_em_ess_fraction
                        ))
  )
}

#' Performs the Data-Adaptive Initialization Strategy via Candidate Model Competition.
#'
#' @param X,y Model matrix and response vector.
#' @param omega_init_guess,sigmaSq_init_guess Initial guesses for omega and sigmaSq.
#'   In the E-step, beta_j | nu_j, lambda_j uses variance omega * nu_j / lambda_j.
#' @param iter_pre_opt Iterations for the pre-optimization EM stage for each candidate.
#' @param pre_opt_burnin,pre_opt_samples Burn-in and saved Gibbs samples within each EM E-step.
#' @param selection_samples Number of posterior samples per candidate used by
#'   the reverse-logistic selection stage.
#' @param IS_period Pre-optimization EM resampling period. 1 runs a fresh Gibbs
#'   E-step every EM iteration; values above 1 reuse samples with importance
#'   weights between resampling iterations.
#' @param min_em_ess_fraction Minimum effective sample size fraction required
#'   to reuse pre-optimization EM samples with importance weights.
#' @param candidates Candidate models and their initial a,b values.
#' @param woodbury Logical, use Woodbury identity in beta updates.
#' @param delta1,delta2,delta3 Convergence tolerances.
#' @param window_size Size of the convergence window.
#' @param diagX Logical, assume diagonal X.
#' @param verbose Logical, print EM and selection progress messages.
#' @return A list containing winning_params and winner_name.
#' @export
initialize_adaptive <- function(X, y,
                                omega_init_guess,
                                sigmaSq_init_guess,
                                iter_pre_opt,
                                candidates, woodbury,
                                pre_opt_burnin = 500,
                                pre_opt_samples = 1000,
                                selection_samples = 200,
                                IS_period = 1,
                                min_em_ess_fraction = 0.1,
                                delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                                window_size = 5, diagX = FALSE,
                                verbose = TRUE) {
  selection_samples <- as.integer(selection_samples)
  if (!is.finite(selection_samples) || selection_samples < 1) {
    stop("selection_samples must be a positive integer.")
  }
  IS_period <- as.integer(IS_period)
  if (!is.finite(IS_period) || IS_period < 1) {
    stop("IS_period must be a positive integer.")
  }
  if (isTRUE(verbose)) {
    cat("Stage 1: Pre-optimizing each candidate model via Gibbs-within-EM...\n")
  }

  pre_optimized_params <- list()
  pre_optimized_states <- list()
  pre_optimized_trajectories <- list()
  engine_args <- list(
    X = X, y = y,
    max_iter = iter_pre_opt,
    iter_burnin = pre_opt_burnin,
    iter_samples = pre_opt_samples,
    delta1 = delta1, delta2 = delta2, delta3 = delta3,
    window_size = window_size, woodbury = woodbury,
    IS_period = IS_period,
    min_em_ess_fraction = min_em_ess_fraction,
    diagX = diagX
  )

  for (name in names(candidates)) {
    if (isTRUE(verbose)) {
      cat("Optimizing candidate model:", name, "...\n")
    }
    candidate_ab <- candidates[[name]]
    opt_res <- do.call(eb_run_em_engine, c(
      list(
        a_init = candidate_ab$a,
        b_init = candidate_ab$b,
        omega_init = omega_init_guess,
        sigmaSq_init = sigmaSq_init_guess,
        verbose = verbose
      ),
      engine_args
    ))
    pre_optimized_params[[name]] <- opt_res$params
    pre_optimized_states[[name]] <- opt_res$final_state
    pre_optimized_trajectories[[name]] <- opt_res$trajectories
  }

  if (isTRUE(verbose)) {
    cat("Stage 2: Pre-sampling each candidate via Gibbs sampler...\n")
  }
  presampled_betas <- list()

  for (name in names(candidates)) {
    params <- pre_optimized_params[[name]]
    starting <- pre_optimized_states[[name]]
    presamples <- eb_getsamples_emp(
      num = selection_samples, X = X, y = y,
      a = params$a, b = params$b,
      omega = params$omega, sigmaSq = params$sigmaSq,
      beta0 = starting$beta0, nu0 = starting$nu0,
      lambda0 = starting$lambda0, xi0 = starting$xi0,
      burn = 0,
      woodbury = woodbury, diagX = diagX
    )
    presampled_betas[[name]] <- presamples$beta
  }

  if (isTRUE(verbose)) {
    cat("Stage 3: Starting model selection ...\n")
  }
  selection_scores <- eb_precompute_selection_scores(
    presampled_betas = presampled_betas,
    pre_optimized_params = pre_optimized_params,
    X = X,
    y = y,
    diagX = diagX
  )

  importance_result <- eb_reverse_logistic_selection(
    presampled_betas = presampled_betas,
    pre_optimized_params = pre_optimized_params,
    X = X,
    y = y,
    diagX = diagX,
    kernel_scores = selection_scores
  )
  avg_probs <- importance_result$model_probabilities
  winner_name <- names(which.max(avg_probs))

  list(
    winning_params = pre_optimized_params[[winner_name]],
    winner_name = winner_name,
    model_probabilities = avg_probs,
    pre_optimized_params = pre_optimized_params,
    pre_optimized_trajectories = pre_optimized_trajectories,
    selection_samples = selection_samples,
    IS_period = IS_period,
    min_em_ess_fraction = min_em_ess_fraction,
    importance = importance_result
  )
}

eb_initial_values <- function(X, y, epsilon = 1e-6,
                              init_method = c("ridge", "olasso"),
                              ridge_lambda = NULL,
                              diagX = FALSE) {
  init_method <- match.arg(init_method)
  n <- nrow(X)

  screen <- if (init_method == "ridge") {
    eb_ridge_screen(
      X = X, y = y, epsilon = epsilon,
      ridge_lambda = ridge_lambda, diagX = diagX
    )
  } else {
    eb_olasso_screen(X, y, epsilon = epsilon)
  }
  sigmaSq_hat <- max(epsilon, screen$sigmaSq)

  y_sq_sum <- sum(y^2)
  trace_XX <- sum(X^2)
  numerator <- y_sq_sum - (n * sigmaSq_hat)
  omega_hat <- if (trace_XX <= 0 || numerator <= 0) {
    epsilon
  } else {
    max(epsilon, numerator / trace_XX)
  }
  list(
    sigmaSq = sigmaSq_hat,
    omega = omega_hat,
    source = screen$source,
    ridge_lambda = screen$ridge_lambda
  )
}

#' Run Empirical Bayes Model Competition for TPB Prior Elicitation
#'
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param iter_pre_opt Iterations for the pre-optimization EM stage for each candidate.
#' @param omega_init_guess,sigmaSq_init_guess Optional initial guesses for omega and sigmaSq.
#'   In the E-step, beta_j | nu_j, lambda_j uses variance omega * nu_j / lambda_j.
#' @param init_method Method used to initialize omega and sigmaSq when explicit
#'   initial guesses are not supplied. "ridge" uses a ridge regression residual
#'   variance; "olasso" uses organic lasso.
#' @param ridge_lambda Ridge penalty used when init_method = "ridge". If NULL,
#'   sqrt(n) is used.
#' @param selection_samples Number of posterior samples per candidate used by
#'   the reverse-logistic selection stage.
#' @param IS_period Pre-optimization EM resampling period. 1 runs a fresh Gibbs
#'   E-step every EM iteration; values above 1 reuse samples with importance
#'   weights between resampling iterations.
#' @param min_em_ess_fraction Minimum effective sample size fraction required
#'   to reuse pre-optimization EM samples with importance weights.
#' @param woodbury Logical, use Woodbury identity in beta updates.
#' @param candidates Candidate models and their initial a,b values.
#' @param pre_opt_burnin,pre_opt_samples Burn-in and saved Gibbs samples within each EM E-step.
#' @param delta1,delta2,delta3 Convergence tolerances.
#' @param window_size Size of the convergence window.
#' @param diagX Logical, assume diagonal X.
#' @param verbose Logical, print EM and selection progress messages.
#' @return A list with winner, winner_name, and raw adaptive competition output.
#' @export
run_model_competition <- function(X, y,
                                  iter_pre_opt = 100,
                                  omega_init_guess = NULL,
                                  sigmaSq_init_guess = NULL,
                                  init_method = c("ridge", "olasso"),
                                  ridge_lambda = NULL,
                                  selection_samples = 200,
                                  IS_period = 1,
                                  min_em_ess_fraction = 0.1,
                                  woodbury = TRUE,
                                  candidates = tpb_default_candidates(),
                                  pre_opt_burnin = 200,
                                  pre_opt_samples = 200,
                                  delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                                  window_size = 5,
                                  diagX = FALSE,
                                  verbose = TRUE) {
  init_method <- match.arg(init_method)
  if (is.null(omega_init_guess) || is.null(sigmaSq_init_guess)) {
    initial_values <- eb_initial_values(
      X = X, y = y,
      init_method = init_method,
      ridge_lambda = ridge_lambda,
      diagX = diagX
    )
    if (is.null(omega_init_guess)) {
      omega_init_guess <- initial_values$omega
    }
    if (is.null(sigmaSq_init_guess)) {
      sigmaSq_init_guess <- initial_values$sigmaSq
    }
  }

  adaptive_result <- initialize_adaptive(
    X = X,
    y = y,
    omega_init_guess = omega_init_guess,
    sigmaSq_init_guess = sigmaSq_init_guess,
    iter_pre_opt = iter_pre_opt,
    selection_samples = selection_samples,
    IS_period = IS_period,
    min_em_ess_fraction = min_em_ess_fraction,
    woodbury = woodbury,
    candidates = candidates,
    pre_opt_burnin = pre_opt_burnin,
    pre_opt_samples = pre_opt_samples,
    delta1 = delta1, delta2 = delta2, delta3 = delta3,
    window_size = window_size,
    diagX = diagX,
    verbose = verbose
  )

  list(
    winner = adaptive_result$winning_params,
    winner_name = adaptive_result$winner_name,
    raw = adaptive_result
  )
}
