#' Run Empirical Bayes Model Competition for TPB Prior Elicitation
#'
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param iter_pre_opt Iterations for the pre-optimization EM stage for each candidate.
#' @param omega_init_guess,sigmaSq_init_guess Optional initial guesses for omega and sigmaSq.
#' @param init_option Initialization method for missing omega/sigmaSq guesses, either "ridge" or "olasso".
#' @param pre_opt_burnin,pre_opt_samples burn-in and sample size for each pre-optimization EM step.
#' @param iter_selection Number of post-burn-in samples for model selection.
#' @param woodbury Logical, use Woodbury identity in beta updates.
#' @param ... Additional arguments passed to initialize_adaptive().
#' @return A list with winner, winner_name, and raw adaptive competition output.
#' @export
run_model_competition <- function(X, y,
                                  iter_pre_opt = 100,
                                  omega_init_guess = NULL,
                                  sigmaSq_init_guess = NULL,
                                  init_option = "ridge",
                                  pre_opt_burnin = 1000,
                                  pre_opt_samples = 1000,
                                  iter_selection = 5000,
                                  woodbury = TRUE,
                                  ...) {
  if (is.null(omega_init_guess) || is.null(sigmaSq_init_guess)) {
    initial_values <- eb_initial_values(
      X = X,
      y = y,
      option = init_option,
      sigmaSq_hat = sigmaSq_init_guess
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
    pre_opt_burnin = pre_opt_burnin,
    pre_opt_samples = pre_opt_samples,
    iter_selection = iter_selection,
    woodbury = woodbury,
    ...
  )

  list(
    winner = adaptive_result$winning_params,
    winner_name = adaptive_result$winner_name,
    raw = adaptive_result
  )
}

#' Performs the Data-Adaptive Initialization Strategy via Candidate Model Competition.
#'
#' @param X,y Model matrix and response vector.
#' @param omega_init_guess,sigmaSq_init_guess Initial guesses for omega and sigmaSq.
#' @param iter_pre_opt Iterations for the pre-optimization EM stage for each candidate.
#' @param pre_opt_burnin,pre_opt_samples burn-in and sample size for each pre-optimization EM step.
#' @param iter_burnin_selection Number of burn-in iterations for the model selection chain.
#' @param iter_selection Number of post-burn-in samples for model selection.
#' @param candidates Candidate models and their fixed a,b values.
#' @param woodbury Logical, use Woodbury identity in beta updates.
#' @param delta1,delta2,delta3 Convergence tolerances.
#' @param window_size Size of the convergence window.
#' @param converge_rule Convergence rule, currently "hyper" or "rho".
#' @param approx Logical placeholder kept for API compatibility.
#' @param option "augmented" or "naive" nu update.
#' @param diagX Logical, assume diagonal X.
#' @return A list containing winning_params and winner_name.
#' @export
initialize_adaptive <- function(X, y,
                                omega_init_guess,
                                sigmaSq_init_guess,
                                iter_pre_opt = 100,
                                pre_opt_burnin = 1000,
                                pre_opt_samples = 1000,
                                iter_burnin_selection = 0,
                                iter_selection = 5000,
                                candidates = list(
                                  horseshoe = list(a = 0.5, b = 0.5),
                                  normal_gamma = list(a = 0.5, b = 10.0),
                                  studentt = list(a = 10.0, b = 1.0)
                                ),
                                woodbury = TRUE,
                                delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                                window_size = 5, converge_rule = "hyper",
                                approx = FALSE, option = "augmented", diagX = FALSE) {
  cat("Stage 1: Pre-optimizing each candidate model...\n")

  pre_optimized_params <- list()
  pre_optimized_states <- list()
  engine_args <- list(
    X = X, y = y,
    max_iter = iter_pre_opt,
    iter_burnin = pre_opt_burnin,
    iter_samples = pre_opt_samples,
    delta1 = delta1, delta2 = delta2, delta3 = delta3,
    window_size = window_size, converge_rule = converge_rule,
    woodbury = woodbury, approx = approx, option = option, diagX = diagX
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
        null.a = FALSE,
        null.b = FALSE,
        null.omega = TRUE,
        null.sigmaSq = TRUE
      ),
      engine_args
    ))
    pre_optimized_params[[name]] <- opt_res$params
    pre_optimized_states[[name]] <- opt_res$final_state
  }

  cat("Stage 2: pre-sampling ...\n")
  presampled_betas <- list()
  presampled_nus <- list()
  presampled_lambdas <- list()

  for (name in names(candidates)) {
    params <- pre_optimized_params[[name]]
    starting <- pre_optimized_states[[name]]
    presamples <- eb_getsamples(
      num = iter_selection, X = X, y = y,
      a = params$a, b = params$b,
      omega = params$omega, sigmaSq = params$sigmaSq,
      beta0 = starting$beta0, nu0 = starting$nu0,
      lambda0 = starting$lambda0, burn = 0,
      xi0 = starting$xi0, woodbury = woodbury,
      approx = approx, diagX = diagX
    )
    presampled_betas[[name]] <- presamples$beta
    presampled_nus[[name]] <- presamples$nu
    presampled_lambdas[[name]] <- presamples$lambda
  }

  cat("Stage 3: Starting iterative model selection ...\n")
  model_names <- names(candidates)
  prob_matrix <- matrix(NA, nrow = iter_selection, ncol = length(candidates))
  current_model_name <- sample(model_names, 1)
  current_beta <- presampled_betas[[current_model_name]][1, ]
  current_nu <- presampled_nus[[current_model_name]][1, ]
  current_lambda <- presampled_lambdas[[current_model_name]][1, ]
  iter_total <- iter_burnin_selection + iter_selection

  for (j in 1:iter_total) {
    log_weights <- sapply(model_names, function(name) {
      eb_calculate_marginal_loglik_beta(
        beta_vec = current_beta,
        model_params = pre_optimized_params[[name]]
      )
    })

    weights <- exp(log_weights - max(log_weights, na.rm = TRUE))
    probs <- weights / sum(weights, na.rm = TRUE)
    current_model_name <- sample(model_names, 1, prob = probs)
    sample_row_idx <- sample(1:iter_selection, 1)
    current_beta <- presampled_betas[[current_model_name]][sample_row_idx, ]
    current_nu <- presampled_nus[[current_model_name]][sample_row_idx, ]
    current_lambda <- presampled_lambdas[[current_model_name]][sample_row_idx, ]

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
    model_probabilities = avg_probs
  )
}

eb_initial_values <- function(X, y, option = c("ridge", "olasso"),
                              sigmaSq_hat = NULL, epsilon = 1e-6) {
  option <- match.arg(option)
  if (is.null(sigmaSq_hat)) {
    sigmaSq_hat <- eb_estimate_sigmaSq(X, y, option = option, epsilon = epsilon)
  }
  sigmaSq_hat <- eb_positive_scalar(sigmaSq_hat, epsilon = epsilon,
                                    name = "sigmaSq")
  omega_hat <- eb_estimate_omega_mom(X, y, sigmaSq_hat, epsilon = epsilon)

  list(omega = omega_hat, sigmaSq = sigmaSq_hat)
}

eb_estimate_sigmaSq <- function(X, y, option = c("ridge", "olasso"),
                                epsilon = 1e-6) {
  option <- match.arg(option)
  y <- as.vector(y)

  if (option == "ridge") {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package 'glmnet' is required for ridge initialization.", call. = FALSE)
    }
    sigmaSq_hat <- eb_estimate_sigmaSq_ridge(X, y, epsilon = epsilon)
  } else {
    if (!requireNamespace("natural", quietly = TRUE)) {
      stop("Package 'natural' is required for olasso initialization.", call. = FALSE)
    }
    olasso <- natural::olasso(X, y, intercept = FALSE)
    sigmaSq_hat <- olasso$sig_obj_1
  }

  eb_positive_scalar(sigmaSq_hat, epsilon = epsilon, name = "sigmaSq")
}

eb_estimate_sigmaSq_ridge <- function(X, y, epsilon = 1e-6) {
  cv_fit <- tryCatch(
    glmnet::cv.glmnet(X, y, alpha = 0, intercept = FALSE, standardize = FALSE),
    error = function(e) e
  )

  if (!inherits(cv_fit, "error")) {
    y_hat <- as.vector(stats::predict(cv_fit, s = "lambda.min", newx = X))
    return(eb_positive_scalar(stats::var(y - y_hat), epsilon = epsilon, name = "sigmaSq"))
  }

  warning(
    sprintf(
      "Ridge CV initialization failed (%s). Falling back to ridge path initialization.",
      conditionMessage(cv_fit)
    )
  )

  ridge_path <- glmnet::glmnet(X, y, alpha = 0, intercept = FALSE, standardize = FALSE)
  pred_mat <- stats::predict(ridge_path, newx = X)
  resid_mat <- sweep(pred_mat, 1, y, FUN = "-")
  resid_var <- apply(resid_mat, 2, stats::var)
  best_idx <- which.min(resid_var)

  eb_positive_scalar(resid_var[best_idx], epsilon = epsilon, name = "sigmaSq")
}

eb_estimate_omega_mom <- function(X, y, sigmaSq_hat, epsilon = 1e-6) {
  n <- nrow(X)
  y_sq_sum <- sum(as.vector(y)^2)
  trace_XX <- sum(X^2)

  if (!is.finite(trace_XX) || trace_XX <= 0) {
    warning("Method of Moments estimate for omega has non-positive trace(X'X). Returning epsilon.")
    return(epsilon)
  }

  numerator <- y_sq_sum - n * sigmaSq_hat
  if (!is.finite(numerator) || numerator <= 0) {
    warning("Method of Moments estimate for omega is non-positive. Returning epsilon.")
    return(epsilon)
  }

  eb_positive_scalar(numerator / trace_XX, epsilon = epsilon, name = "omega")
}

eb_positive_scalar <- function(x, epsilon = 1e-6, name = "value") {
  if (length(x) != 1 || !is.finite(x) || x <= 0) {
    warning(sprintf("Initialization estimate for %s is not positive and finite. Returning epsilon.", name))
    return(epsilon)
  }
  max(x, epsilon)
}

eb_run_em_engine <- function(X, y,
                             a_init, b_init, omega_init, sigmaSq_init,
                             null.a, null.b, null.omega, null.sigmaSq,
                             max_iter, iter_burnin, iter_samples,
                             delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                             window_size = 5, converge_rule = "hyper",
                             woodbury = FALSE, approx = FALSE,
                             option = "augmented", diagX = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  a_vec <- c(a_init)
  b_vec <- c(b_init)
  omega_vec <- c(omega_init)
  sigmaSq_vec <- c(sigmaSq_init)
  beta_diff_vec <- numeric(0)
  last_beta_hat <- NULL
  beta0 <- rnorm(p, mean = 0, sd = sqrt(sigmaSq_init))
  nu0 <- rgamma(p, 1, 1)
  lambda0 <- rgamma(p, 1, 1)
  xi0 <- if (option == "augmented") {
    a_init * sqrt(exp(log(beta0^2) + log(lambda0) - log(omega_init)) / 2)
  } else {
    NULL
  }

  for (k in 1:max_iter) {
    samples <- eb_getsamples(
      num = iter_samples, X = X, y = y,
      a = a_vec[k], b = b_vec[k],
      omega = omega_vec[k], sigmaSq = sigmaSq_vec[k],
      beta0 = beta0, nu0 = nu0, lambda0 = lambda0,
      burn = iter_burnin, xi0 = xi0,
      woodbury = woodbury, approx = approx, diagX = diagX
    )

    a_new <- if (null.a) eb_M_step(samples$nu, a_vec[k]) else a_vec[k]
    b_new <- if (null.b) eb_M_step(samples$lambda, b_vec[k]) else b_vec[k]
    sigmaSq_new <- if (null.sigmaSq) eb_M_step_sigmaSq(samples$beta, X, y, n, diagX) else sigmaSq_vec[k]
    omega_new <- if (null.omega) eb_M_step_phi(samples$beta, samples$lambda, samples$nu) else omega_vec[k]

    beta0 <- samples$beta[iter_samples, ]
    nu0 <- samples$nu[iter_samples, ]
    lambda0 <- samples$lambda[iter_samples, ]
    xi0 <- if (!is.null(samples$xi)) samples$xi[iter_samples, ] else NULL
    current_beta_hat <- colMeans(samples$beta)
    beta_diff_vec <- c(beta_diff_vec, if (is.null(last_beta_hat)) NA else mean((current_beta_hat - last_beta_hat)^2))
    last_beta_hat <- current_beta_hat
    a_vec <- c(a_vec, a_new)
    b_vec <- c(b_vec, b_new)
    sigmaSq_vec <- c(sigmaSq_vec, sigmaSq_new)
    omega_vec <- c(omega_vec, omega_new)

    converged <- if (converge_rule == "hyper") {
      eb_isConverged_hyper(a_vec, b_vec, sigmaSq_vec, omega_vec,
                           delta1 = delta1, delta2 = delta2,
                           window_size = window_size)
    } else {
      eb_isConverged_rho(a_vec, b_vec, omega_vec,
                         delta = delta3, window_size = window_size)
    }
    if (converged) {
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
    final_state = list(beta0 = beta0, nu0 = nu0,
                       lambda0 = lambda0, xi0 = xi0),
    trajectories = list(a_traj = a_vec, b_traj = b_vec,
                        sigmaSq_traj = sigmaSq_vec,
                        omega_traj = omega_vec,
                        beta_diff_traj = beta_diff_vec)
  )
}

eb_getsamples <- function(num, X, y, a, b, omega, sigmaSq,
                          beta0, nu0, lambda0, burn, xi0 = NULL,
                          woodbury = FALSE, approx = FALSE,
                          diagX = FALSE) {
  aug <- !is.null(xi0)
  p <- ncol(X)
  beta_samples <- matrix(NA, nrow = num + burn, ncol = p)
  nu_samples <- matrix(NA, nrow = num + burn, ncol = p)
  lambda_samples <- matrix(NA, nrow = num + burn, ncol = p)
  beta_samples[1, ] <- beta0
  nu_samples[1, ] <- nu0
  lambda_samples[1, ] <- lambda0
  if (aug) {
    xi_samples <- matrix(NA, nrow = num + burn, ncol = p)
    xi_samples[1, ] <- xi0
  }

  for (i in 2:(num + burn)) {
    beta_samples[i, ] <- eb_Gibbs_beta(
      X = X, y = y, omega = omega, sigmaSq = sigmaSq,
      nu = nu_samples[i - 1, ], lambda = lambda_samples[i - 1, ],
      woodbury = woodbury, approx = approx, diagX = diagX
    )
    if (aug) {
      nu_samples[i, ] <- eb_aug_Gibbs_nu(
        p = p, omega = omega, beta = beta_samples[i, ],
        lambda = lambda_samples[i - 1, ],
        xi = xi_samples[i - 1, ], a = a
      )
      xi_samples[i, ] <- eb_Gibbs_xi(a = a, p = p, nu = nu_samples[i, ])
    } else {
      nu_samples[i, ] <- eb_Naive_Gibbs_nu(
        a = a, p = p, omega = omega,
        beta = beta_samples[i, ],
        lambda = lambda_samples[i - 1, ]
      )
    }
    lambda_samples[i, ] <- eb_Gibbs_lambda(
      b = b, p = p, omega = omega,
      beta = beta_samples[i, ], nu = nu_samples[i, ]
    )
  }

  res <- list(
    beta = beta_samples[(burn + 1):(burn + num), ],
    nu = nu_samples[(burn + 1):(burn + num), ],
    lambda = lambda_samples[(burn + 1):(burn + num), ]
  )
  if (aug) {
    res$xi <- xi_samples[(burn + 1):(burn + num), ]
  }
  res
}

eb_Gibbs_beta <- function(X, y, omega, sigmaSq, nu, lambda,
                          woodbury = FALSE, approx = FALSE,
                          diagX = FALSE) {
  if (!diagX) {
    if (!woodbury) {
      eta <- omega * nu / lambda
      Xrteta <- (X / sqrt(sigmaSq)) %*% diag(sqrt(eta))
      AA <- eigen(t(Xrteta) %*% Xrteta)
      nonnegative_eigenVals <- ifelse(AA$values + 1 > 0, 1 / (AA$values + 1), 0)
      inverse <- diag(sqrt(eta)) %*% AA$vectors %*%
        diag(nonnegative_eigenVals) %*% t(AA$vectors) %*% diag(sqrt(eta))
      mean_beta <- inverse %*% t(X) %*% (y / sigmaSq)
      beta <- diag(sqrt(eta)) %*% AA$vectors %*%
        diag(sqrt(nonnegative_eigenVals)) %*% rnorm(length(mean_beta)) + mean_beta
    } else {
      n <- nrow(X)
      p <- ncol(X)
      d <- omega * nu / lambda
      u <- rnorm(p, 0, sqrt(d))
      delta <- rnorm(n, 0, 1)
      XD <- sweep(X, 2, d, FUN = "*")
      M <- XD %*% t(X) + sigmaSq * diag(n)
      eig <- eigen(M, symmetric = TRUE)
      M_eig_vals_inv <- ifelse(eig$values > 1e-8, 1 / eig$values, 0)
      M_inv <- eig$vectors %*% diag(M_eig_vals_inv) %*% t(eig$vectors)
      rhs <- y - X %*% u - sqrt(sigmaSq) * delta
      w <- M_inv %*% rhs
      beta <- u + sweep(t(XD) %*% w, 1, 1, FUN = "*")
    }
    return(beta)
  }

  mean_beta <- ((diag(X) * y) / sigmaSq) /
    (diag(X)^2 / sigmaSq + lambda / (omega * nu))
  var_beta <- 1 / (diag(X)^2 / sigmaSq + lambda / (omega * nu))
  rnorm(length(mean_beta)) * sqrt(var_beta) + mean_beta
}

eb_aug_Gibbs_nu <- function(p, omega, beta, lambda, xi, a) {
  chi <- exp(log(beta^2) + log(lambda) - log(omega))
  statmod::rinvgauss(p, mean = sqrt((chi + 2 * xi) / (2 * a)),
                     shape = chi + 2 * xi)
}

eb_Naive_Gibbs_nu <- function(a, p, omega, beta, lambda) {
  vapply(1:p, function(d) {
    chi <- max(exp(log(beta[d]^2) + log(lambda[d]) - log(omega)), 1e-306)
    GIGrvg::rgig(1, lambda = a - 0.5, chi = chi, psi = 2 * a)
  }, numeric(1))
}

eb_Gibbs_xi <- function(a, p, nu) {
  rgamma(p, shape = a, rate = 1 / nu)
}

eb_Gibbs_lambda <- function(b, p, omega, beta, nu) {
  rgamma(p, shape = b + 0.5, rate = b + 0.5 * beta^2 / (omega * nu))
}

eb_M_step <- function(sample_matrix, param_current) {
  empMean_log <- colMeans(log(sample_matrix + .Machine$double.eps))
  idigamma_input <- log(param_current) + mean(empMean_log)
  if (!is.finite(idigamma_input)) {
    warning("M-step for a/b received non-finite input.")
    return(1e-4)
  }
  max(1e-6, eb_idigamma(idigamma_input))
}

eb_idigamma <- function(y) {
  lower <- 1e-8
  upper <- 1
  while (digamma(upper) < y && upper < 1e8) {
    upper <- upper * 2
  }
  stats::uniroot(function(x) digamma(x) - y,
                 lower = lower, upper = upper)$root
}

eb_M_step_sigmaSq <- function(beta_matrix, X, y, n, diagX = FALSE) {
  num_samples <- nrow(beta_matrix)
  if (!diagX) {
    ssr_per_sample <- vapply(1:num_samples, function(i) {
      residuals_i <- y - X %*% beta_matrix[i, ]
      sum(residuals_i^2)
    }, numeric(1))
  } else {
    d <- diag(X)
    ssr_per_sample <- vapply(1:num_samples, function(i) {
      residuals_i <- y - d * beta_matrix[i, ]
      sum(residuals_i^2)
    }, numeric(1))
  }
  max(mean(ssr_per_sample) / n, 1e-16)
}

eb_M_step_phi <- function(beta_matrix, lambda_matrix, nu_matrix) {
  empMean <- colMeans(beta_matrix^2 * lambda_matrix /
                        (nu_matrix + .Machine$double.eps))
  if (any(!is.finite(empMean))) {
    warning("Non-finite values detected in omega M-step. Using finite values only.")
    empMean <- empMean[is.finite(empMean)]
    if (length(empMean) == 0) {
      return(1e-4)
    }
  }
  mean(empMean)
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

eb_isConverged_rho <- function(a_vec, b_vec, omega_vec,
                               delta, window_size = 3) {
  if (length(a_vec) < window_size + 1) {
    return(FALSE)
  }
  max(abs(diff(tail(a_vec, window_size + 1))),
      abs(diff(tail(b_vec, window_size + 1))),
      abs(diff(tail(omega_vec, window_size + 1))),
      na.rm = TRUE) < delta
}

eb_calculate_marginal_loglik_beta <- function(beta_vec, model_params) {
  a <- model_params$a
  b <- model_params$b
  phi <- model_params$omega * b / a
  loglik_individual <- log(pmax(
    eb_d_tpb(beta_vec = beta_vec, a = a, b = b, phi = phi),
    .Machine$double.xmin
  ))
  if (any(!is.finite(loglik_individual))) {
    return(-Inf)
  }
  sum(loglik_individual)
}

eb_d_tpb <- function(beta_vec, a, b, phi) {
  vapply(beta_vec, function(beta_i) {
    const <- gamma(0.5 + b) * gamma(a + b) /
      (gamma(a) * gamma(b) * sqrt(2 * pi * phi))
    z_val <- beta_i^2 / (2 * phi)
    U_val <- eb_tricomi_U(0.5 + b, 1.5 - a, z_val)
    if (is.na(U_val)) {
      return(NA_real_)
    }
    const * U_val
  }, numeric(1))
}

eb_tricomi_U <- function(a, b, z) {
  U_val <- tryCatch({
    gsl::hyperg_U(a, b, z)
  }, error = function(e) NA_real_)
  return(U_val)
}
