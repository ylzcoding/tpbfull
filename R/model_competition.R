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

eb_olasso_active_prop <- function(ol_fit, epsilon = 1e-8) {
  if (is.null(ol_fit) || is.null(ol_fit$beta_1)) {
    return(NA_real_)
  }
  mean(abs(ol_fit$beta_1) > epsilon)
}

eb_olasso_active_count <- function(ol_fit, epsilon = 1e-8) {
  if (is.null(ol_fit) || is.null(ol_fit$beta_1)) {
    return(NA_integer_)
  }
  sum(abs(ol_fit$beta_1) > epsilon)
}

eb_candidate_model_log_prior <- function(X, y, candidates,
                                         active_prop = NULL,
                                         active_count = NULL,
                                         sparsity_threshold = 0.10,
                                         model_prior_method = c("complexity", "binomial"),
                                         model_size_prior_alpha = 1,
                                         model_size_prior_beta = 1) {
  model_prior_method <- match.arg(model_prior_method)
  if (is.null(active_prop) || !is.finite(active_prop) ||
      is.null(active_count) || !is.finite(active_count)) {
    ol_summary <- tryCatch({
      ol_fit <- natural::olasso(X, y, intercept = FALSE)
      list(
        active_prop = eb_olasso_active_prop(ol_fit),
        active_count = eb_olasso_active_count(ol_fit)
      )
    }, error = function(e) list(active_prop = NA_real_, active_count = NA_real_))
    if (is.null(active_prop) || !is.finite(active_prop)) {
      active_prop <- ol_summary$active_prop
    }
    if (is.null(active_count) || !is.finite(active_count)) {
      active_count <- ol_summary$active_count
    }
  }
  log_prior <- stats::setNames(rep(0, length(candidates)), names(candidates))
  if (!is.finite(active_prop) || !is.finite(active_count)) {
    return(list(log_prior = log_prior, active_prop = active_prop,
                active_count = active_count,
                source = "olasso",
                sparsity_threshold = sparsity_threshold,
                model_prior_method = model_prior_method,
                model_size_prior_alpha = model_size_prior_alpha,
                model_size_prior_beta = model_size_prior_beta))
  }

  candidate_ab <- do.call(rbind, lapply(candidates, function(x) {
    c(a = x$a, b = x$b)
  }))
  candidate_names <- rownames(candidate_ab)
  is_hs <- candidate_names == "hs" |
    (candidate_ab[, "a"] <= 0.5 & candidate_ab[, "b"] <= 0.5)
  is_dense <- candidate_names %in% c("normal_gamma", "ridge", "student_t") |
    candidate_ab[, "b"] >= 5 |
    candidate_ab[, "a"] >= 5 |
    (candidate_ab[, "a"] >= 5 & candidate_ab[, "b"] >= 5)

  p <- ncol(X)
  sparse_cutoff <- floor(sparsity_threshold * p)

  if (model_prior_method == "complexity") {
    if (active_count <= sparse_cutoff) {
      dense_gap <- sparse_cutoff + 1 - active_count
      log_prior[is_dense] <- -dense_gap * log(p)
    } else {
      hs_gap <- active_count - sparse_cutoff
      log_prior[is_hs] <- -hs_gap * log(p)
    }
    if (any(is_hs)) {
      log_prior[is_hs] <- log_prior[is_hs] - log(sum(is_hs))
    }
    if (any(is_dense)) {
      log_prior[is_dense] <- log_prior[is_dense] - log(sum(is_dense))
    }
  } else {
    q_hat <- (active_count + model_size_prior_alpha) /
      (p + model_size_prior_alpha + model_size_prior_beta)
    q_hat <- min(max(q_hat, .Machine$double.eps), 1 - .Machine$double.eps)

    log_p_sparse <- stats::pbinom(
      sparse_cutoff, size = p, prob = q_hat, log.p = TRUE
    )
    log_p_dense <- stats::pbinom(
      sparse_cutoff, size = p, prob = q_hat,
      lower.tail = FALSE, log.p = TRUE
    )

    if (any(is_hs)) {
      log_prior[is_hs] <- log_p_sparse - log(sum(is_hs))
    }
    if (any(is_dense)) {
      log_prior[is_dense] <- log_p_dense - log(sum(is_dense))
    }
  }

  list(log_prior = log_prior, active_prop = active_prop,
       active_count = active_count,
       source = "olasso", sparsity_threshold = sparsity_threshold,
       model_prior_method = model_prior_method,
       model_size_prior_alpha = model_size_prior_alpha,
       model_size_prior_beta = model_size_prior_beta)
}

eb_calculate_marginal_loglik_beta <- function(beta_vec, model_params) {
  a <- model_params$a
  b <- model_params$b
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

  vapply(beta_vec, function(beta_i) {
    z_val <- beta_i^2 / (2 * phi)
    if (!is.finite(z_val) || z_val < 0) {
      return(NA_real_)
    }
    U_val <- tryCatch({
      gsl::hyperg_U(U_a, U_b, z_val)
    }, error = function(e) NA_real_)

    if (!is.finite(U_val) || U_val <= 0) {
      U_val <- tricomi_U_integral(U_a, U_b, z_val)
    }

    if (!is.finite(U_val) || U_val <= 0) {
      return(NA_real_)
    }
    exp(pmin(log_const + log(U_val), log(.Machine$double.xmax)))
  }, numeric(1))
}

eb_getsamples_emp <- function(num, X, y, a, b, omega, sigmaSq,
                              beta0, psi0, zeta0, burn = 0,
                              woodbury = FALSE, diagX = FALSE) {
  total <- num + burn
  if (total < 1) {
    stop("num + burn must be at least 1.")
  }

  p <- ncol(X)
  beta_samples <- matrix(NA_real_, nrow = total, ncol = p)
  psi_samples <- matrix(NA_real_, nrow = total, ncol = p)
  zeta_samples <- matrix(NA_real_, nrow = total, ncol = p)

  beta_samples[1, ] <- beta0
  psi_samples[1, ] <- pmax(psi0, 1e-16)
  zeta_samples[1, ] <- pmax(zeta0, 1e-16)

  if (total >= 2) {
    for (i in 2:total) {
      phi_current <- eb_positive_finite(omega * b / a, 1e-16)
      beta_samples[i, ] <- Gibbs_beta(
        X = X, y = y, a = a, b = b,
        phi = phi_current,
        sigmaSq = eb_positive_finite(sigmaSq, 1e-16),
        psi = psi_samples[i - 1, ],
        woodbury = woodbury, diagX = diagX
      )
      zeta_samples[i, ] <- Gibbs_zeta(
        p = p, a = a, b = b, psi = psi_samples[i - 1, ]
      )
      psi_samples[i, ] <- Gibbs_psi(
        p = p, b = b, phi = phi_current,
        beta = beta_samples[i, ], zeta = zeta_samples[i, ]
      )
    }
  }

  keep <- (burn + 1):total
  list(
    beta = beta_samples[keep, , drop = FALSE],
    psi = psi_samples[keep, , drop = FALSE],
    zeta = zeta_samples[keep, , drop = FALSE]
  )
}

# Current Gibbs functions use the IG-IG TPB representation:
#   beta | psi, phi ~ N(0, phi * psi), phi = omega * b / a.
# This is equivalent to the EM formulas with
#   a * nu = 1 / zeta and b * lambda = 1 / (zeta * psi).
eb_m_step_a <- function(zeta_matrix) {
  target <- -mean(log(pmax(zeta_matrix, .Machine$double.xmin)))
  candidate <- tryCatch(bzinb::idigamma(target), error = function(e) NA_real_)
  eb_positive_finite(as.numeric(candidate), 1e-6)
}

eb_m_step_b <- function(psi_matrix, zeta_matrix) {
  target <- -mean(log(pmax(psi_matrix, .Machine$double.xmin))) -
    mean(log(pmax(zeta_matrix, .Machine$double.xmin)))
  candidate <- tryCatch(bzinb::idigamma(target), error = function(e) NA_real_)
  eb_positive_finite(as.numeric(candidate), 1e-6)
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

eb_m_step_omega <- function(beta_matrix, psi_matrix, a_current, b_current) {
  vals <- beta_matrix^2 * a_current / (b_current * pmax(psi_matrix, 1e-16))
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
                             update_a = TRUE,
                             update_b = TRUE,
                             update_omega = TRUE,
                             update_sigmaSq = TRUE,
                             diagX = FALSE) {
  p <- ncol(X)

  a_vec <- c(eb_positive_finite(a_init, 1e-6))
  b_vec <- c(eb_positive_finite(b_init, 1e-6))
  omega_vec <- c(eb_positive_finite(omega_init, 1e-16))
  sigmaSq_vec <- c(eb_positive_finite(sigmaSq_init, 1e-16))
  beta_diff_vec <- numeric(0)

  beta0 <- rnorm(p, mean = 0, sd = sqrt(sigmaSq_vec[1]))
  nu0 <- rgamma(p, 1, 1)
  lambda0 <- rgamma(p, 1, 1)
  zeta0 <- 1 / pmax(a_vec[1] * nu0, 1e-16)
  psi0 <- b_vec[1] * nu0 / pmax(a_vec[1] * lambda0, 1e-16)
  last_beta_hat <- NULL

  for (k in seq_len(max_iter)) {
    samples <- eb_getsamples_emp(
      num = iter_samples, X = X, y = y,
      a = a_vec[k], b = b_vec[k],
      omega = omega_vec[k], sigmaSq = sigmaSq_vec[k],
      beta0 = beta0, psi0 = psi0, zeta0 = zeta0,
      burn = iter_burnin,
      woodbury = woodbury, diagX = diagX
    )

    a_new <- if (update_a) eb_m_step_a(samples$zeta) else a_vec[k]
    b_new <- if (update_b) eb_m_step_b(samples$psi, samples$zeta) else b_vec[k]
    sigmaSq_new <- if (update_sigmaSq) eb_m_step_sigmaSq(samples$beta, X, y, diagX) else sigmaSq_vec[k]
    omega_new <- if (update_omega) {
      eb_m_step_omega(samples$beta, samples$psi, a_vec[k], b_vec[k])
    } else {
      omega_vec[k]
    }

    beta0 <- samples$beta[nrow(samples$beta), ]
    psi0 <- samples$psi[nrow(samples$psi), ]
    zeta0 <- samples$zeta[nrow(samples$zeta), ]

    beta_hat <- colMeans(samples$beta)
    beta_diff_vec <- c(
      beta_diff_vec,
      if (is.null(last_beta_hat)) NA_real_ else mean((beta_hat - last_beta_hat)^2)
    )
    last_beta_hat <- beta_hat

    a_vec <- c(a_vec, a_new)
    b_vec <- c(b_vec, b_new)
    sigmaSq_vec <- c(sigmaSq_vec, sigmaSq_new)
    omega_vec <- c(omega_vec, omega_new)

    hyper_converged <- eb_isConverged_hyper(
      a_vec, b_vec, sigmaSq_vec, omega_vec,
      delta1 = delta1, delta2 = delta2, window_size = window_size
    )
    beta_converged <- length(beta_diff_vec) >= window_size &&
      max(tail(beta_diff_vec, window_size), na.rm = TRUE) < delta3

    if (hyper_converged && beta_converged) {
      cat("Engine converged after", k, "iterations.\n")
      break
    }

    cat(k, "-th iteration completed.\n")
    cat("a=", a_vec[k + 1], "\n")
    cat("b=", b_vec[k + 1], "\n")
    cat("omega=", omega_vec[k + 1], "\n")
    cat("sigmaSq=", sigmaSq_vec[k + 1], "\n")
  }

  l <- length(a_vec)
  list(
    params = list(a = a_vec[l], b = b_vec[l],
                  omega = omega_vec[l], sigmaSq = sigmaSq_vec[l]),
    final_state = list(beta0 = beta0, psi0 = psi0, zeta0 = zeta0),
    trajectories = list(a_traj = a_vec, b_traj = b_vec,
                        sigmaSq_traj = sigmaSq_vec,
                        omega_traj = omega_vec,
                        beta_diff_traj = beta_diff_vec)
  )
}

#' Performs the Data-Adaptive Initialization Strategy via Candidate Model Competition.
#'
#' @param X,y Model matrix and response vector.
#' @param omega_init_guess,sigmaSq_init_guess Initial guesses for omega and sigmaSq.
#'   The Gibbs E-step uses phi = omega * b / a.
#' @param iter_pre_opt Iterations for the pre-optimization EM stage for each candidate.
#' @param pre_opt_burnin,pre_opt_samples Burn-in and saved Gibbs samples within each EM E-step.
#' @param iter_burnin_selection Burn-in iterations for the model indicator chain.
#' @param iter_selection Number of post-burn-in samples for model selection.
#' @param candidates Candidate models and their initial a,b values.
#' @param woodbury Logical, use Woodbury identity in beta updates.
#' @param update_a,update_b Logical, update TPB shape parameters during pre-optimization.
#' @param update_omega,update_sigmaSq Logical, update omega and sigmaSq during pre-optimization.
#' @param delta1,delta2,delta3 Convergence tolerances.
#' @param window_size Size of the convergence window.
#' @param diagX Logical, assume diagonal X.
#' @return A list containing winning_params and winner_name.
#' @export
initialize_adaptive <- function(X, y,
                                omega_init_guess,
                                sigmaSq_init_guess,
                                iter_pre_opt,
                                iter_selection,
                                candidates, woodbury,
                                model_log_prior = NULL,
                                pre_opt_burnin = 500,
                                pre_opt_samples = 1000,
                                iter_burnin_selection = 0,
                                update_a = FALSE,
                                update_b = FALSE,
                                update_omega = TRUE,
                                update_sigmaSq = TRUE,
                                delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                                window_size = 5, diagX = FALSE) {
  cat("Stage 1: Pre-optimizing each candidate model via Gibbs-within-EM...\n")

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
    diagX = diagX
  )

  for (name in names(candidates)) {
    cat("Optimizing candidate model:", name, "...\n")
    candidate_ab <- candidates[[name]]
    opt_res <- do.call(eb_run_em_engine, c(
      list(
        a_init = candidate_ab$a,
        b_init = candidate_ab$b,
        omega_init = omega_init_guess,
        sigmaSq_init = sigmaSq_init_guess,
        update_a = update_a,
        update_b = update_b,
        update_omega = update_omega,
        update_sigmaSq = update_sigmaSq
      ),
      engine_args
    ))
    pre_optimized_params[[name]] <- opt_res$params
    pre_optimized_states[[name]] <- opt_res$final_state
    pre_optimized_trajectories[[name]] <- opt_res$trajectories
  }

  cat("Stage 2: Pre-sampling each candidate via Gibbs sampler...\n")
  presampled_betas <- list()

  for (name in names(candidates)) {
    params <- pre_optimized_params[[name]]
    starting <- pre_optimized_states[[name]]
    presamples <- eb_getsamples_emp(
      num = iter_selection, X = X, y = y,
      a = params$a, b = params$b,
      omega = params$omega, sigmaSq = params$sigmaSq,
      beta0 = starting$beta0, psi0 = starting$psi0,
      zeta0 = starting$zeta0, burn = 0,
      woodbury = woodbury, diagX = diagX
    )
    presampled_betas[[name]] <- presamples$beta
  }

  cat("Stage 3: Starting iterative model selection ...\n")
  model_names <- names(candidates)
  if (is.null(model_log_prior)) {
    model_log_prior <- stats::setNames(rep(0, length(model_names)), model_names)
  } else {
    model_log_prior <- model_log_prior[model_names]
    model_log_prior[!is.finite(model_log_prior)] <- 0
  }
  prob_matrix <- matrix(NA_real_, nrow = iter_selection, ncol = length(candidates))
  colnames(prob_matrix) <- model_names
  current_model_name <- sample(model_names, 1)
  current_beta <- presampled_betas[[current_model_name]][1, ]

  iter_total <- iter_burnin_selection + iter_selection
  for (j in seq_len(iter_total)) {
    log_weights <- sapply(model_names, function(name) {
      eb_calculate_marginal_loglik_beta(
        beta_vec = current_beta,
        model_params = pre_optimized_params[[name]]
      ) + model_log_prior[[name]]
    })

    max_log <- max(log_weights, na.rm = TRUE)
    if (!is.finite(max_log)) {
      probs <- rep(1 / length(model_names), length(model_names))
    } else {
      weights <- exp(log_weights - max_log)
      weight_sum <- sum(weights, na.rm = TRUE)
      probs <- if (is.finite(weight_sum) && weight_sum > 0) {
        weights / weight_sum
      } else {
        rep(1 / length(model_names), length(model_names))
      }
    }

    current_model_name <- sample(model_names, 1, prob = probs)
    sample_row_idx <- sample.int(iter_selection, 1)
    current_beta <- presampled_betas[[current_model_name]][sample_row_idx, ]

    if (j > iter_burnin_selection) {
      prob_matrix[j - iter_burnin_selection, ] <- probs
    }
  }

  avg_probs <- colMeans(prob_matrix)
  names(avg_probs) <- model_names
  winner_name <- names(which.max(avg_probs))

  list(
    winning_params = pre_optimized_params[[winner_name]],
    winner_name = winner_name,
    model_probabilities = avg_probs,
    pre_optimized_params = pre_optimized_params,
    pre_optimized_trajectories = pre_optimized_trajectories
  )
}

eb_initial_values <- function(X, y, epsilon = 1e-6) {
  n <- nrow(X)

  ol_fit <- natural::olasso(X, y, intercept = FALSE)
  sigmaSq_hat <- max(epsilon, ol_fit$sig_obj_1)
  active_prop <- eb_olasso_active_prop(ol_fit)
  active_count <- eb_olasso_active_count(ol_fit)

  y_sq_sum <- sum(y^2)
  trace_XX <- sum(X^2)
  numerator <- y_sq_sum - (n * sigmaSq_hat)
  omega_hat <- if (trace_XX <= 0 || numerator <= 0) {
    epsilon
  } else {
    max(epsilon, numerator / trace_XX)
  }
  list(sigmaSq = sigmaSq_hat, omega = omega_hat,
       active_prop = active_prop, active_count = active_count)
}

#' Run Empirical Bayes Model Competition for TPB Prior Elicitation
#'
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param iter_pre_opt Iterations for the pre-optimization EM stage for each candidate.
#' @param omega_init_guess,sigmaSq_init_guess Optional initial guesses for omega and sigmaSq.
#'   The Gibbs E-step uses phi = omega * b / a.
#' @param iter_selection Number of post-burn-in samples for model selection.
#' @param woodbury Logical, use Woodbury identity in beta updates.
#' @param candidates Candidate models and their initial a,b values.
#' @param pre_opt_burnin,pre_opt_samples Burn-in and saved Gibbs samples within each EM E-step.
#' @param iter_burnin_selection Burn-in iterations for the model indicator chain.
#' @param use_model_prior Logical, add the organic-lasso regime prior to model weights.
#' @param sparsity_threshold Active-proportion threshold for sparse vs dense regimes.
#' @param model_prior_method "complexity" penalizes regime mismatch on a log(p)
#'   model-complexity scale; "binomial" uses the organic-lasso model size to
#'   induce a binomial regime prior.
#' @param model_size_prior_alpha,model_size_prior_beta Beta smoothing constants for
#'   the organic-lasso inclusion probability in the size-based regime prior.
#' @param delta1,delta2,delta3 Convergence tolerances.
#' @param window_size Size of the convergence window.
#' @param diagX Logical, assume diagonal X.
#' @return A list with winner, winner_name, and raw adaptive competition output.
#' @export
run_model_competition <- function(X, y,
                                  iter_pre_opt = 100,
                                  omega_init_guess = NULL,
                                  sigmaSq_init_guess = NULL,
                                  iter_selection = 5000,
                                  woodbury = TRUE,
                                  candidates = list(
                                    hs = list(a = 0.5, b = 0.5),
                                    normal_gamma = list(a = 1, b = 10.0),
                                    ridge = list(a = 10.0, b = 10.0)
                                  ),
                                  pre_opt_burnin = 200,
                                  pre_opt_samples = 200,
                                  iter_burnin_selection = 0,
                                  use_model_prior = TRUE,
                                  sparsity_threshold = 0.10,
                                  model_prior_method = c("complexity", "binomial"),
                                  model_size_prior_alpha = 1,
                                  model_size_prior_beta = 1,
                                  delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                                  window_size = 5,
                                  diagX = FALSE) {
  model_prior_method <- match.arg(model_prior_method)
  initial_values <- NULL
  active_prop <- NA_real_
  active_count <- NA_real_
  if (is.null(omega_init_guess) || is.null(sigmaSq_init_guess) ||
      isTRUE(use_model_prior)) {
    initial_values <- eb_initial_values(X = X, y = y)
    active_prop <- initial_values$active_prop
    active_count <- initial_values$active_count
    if (is.null(omega_init_guess)) {
      omega_init_guess <- initial_values$omega
    }
    if (is.null(sigmaSq_init_guess)) {
      sigmaSq_init_guess <- initial_values$sigmaSq
    }
  }

  prior_info <- if (isTRUE(use_model_prior)) {
    eb_candidate_model_log_prior(
      X = X,
      y = y,
      candidates = candidates,
      active_prop = active_prop,
      active_count = active_count,
      sparsity_threshold = sparsity_threshold,
      model_prior_method = model_prior_method,
      model_size_prior_alpha = model_size_prior_alpha,
      model_size_prior_beta = model_size_prior_beta
    )
  } else {
    list(
      log_prior = stats::setNames(rep(0, length(candidates)), names(candidates)),
      active_prop = NA_real_,
      active_count = NA_real_,
      source = NA_character_,
      sparsity_threshold = sparsity_threshold,
      model_prior_method = model_prior_method,
      model_size_prior_alpha = model_size_prior_alpha,
      model_size_prior_beta = model_size_prior_beta
    )
  }

  adaptive_result <- initialize_adaptive(
    X = X,
    y = y,
    omega_init_guess = omega_init_guess,
    sigmaSq_init_guess = sigmaSq_init_guess,
    iter_pre_opt = iter_pre_opt,
    iter_selection = iter_selection,
    woodbury = woodbury,
    candidates = candidates,
    model_log_prior = prior_info$log_prior,
    pre_opt_burnin = pre_opt_burnin,
    pre_opt_samples = pre_opt_samples,
    iter_burnin_selection = iter_burnin_selection,
    delta1 = delta1, delta2 = delta2, delta3 = delta3,
    window_size = window_size,
    diagX = diagX
  )

  list(
    winner = adaptive_result$winning_params,
    winner_name = adaptive_result$winner_name,
    raw = adaptive_result,
    model_prior = prior_info
  )
}
