#' Run Empirical Bayes Model Competition for TPB Prior Elicitation
#'
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param iter_pre_opt Iterations for the pre-optimization EM stage for each candidate.
#' @param omega_init_guess,sigmaSq_init_guess Optional initial guesses for omega and sigmaSq.
#' @param init_option Initialization method for missing omega/sigmaSq guesses, either "ridge" or "olasso".
#' @param iter_selection Number of post-burn-in samples for model selection.
#' @param woodbury Logical, use Woodbury identity in beta updates.
#' @param candidates Candidate models and their initial a,b values.
#' @param delta1,delta2,delta3 Convergence tolerances.
#' @param window_size Size of the convergence window.
#' @param diagX Logical, assume diagonal X.
#' @return A list with winner, winner_name, and raw adaptive competition output.
#' @export
run_model_competition <- function(X, y,
                                  iter_pre_opt = 100,
                                  omega_init_guess = NULL,
                                  sigmaSq_init_guess = NULL,
                                  init_option = "ridge",
                                  iter_selection = 5000,
                                  woodbury = TRUE,
                                  candidates = list(
                                    horseshoe = list(a = 0.5, b = 0.5),
                                    normal_gamma = list(a = 0.5, b = 5.0),
                                    studentt = list(a = 5.0, b = 0.5)
                                  ),
                                  delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                                  window_size = 5,
                                  diagX = FALSE) {
  if (is.null(omega_init_guess) || is.null(sigmaSq_init_guess)) {
    initial_values <- eb_initial_values(
      X = X,
      y = y,
      option = init_option
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
    iter_selection = iter_selection,
    woodbury = woodbury,
    candidates = candidates,
    delta1 = delta1, delta2 = delta2, delta3 = delta3,
    window_size = window_size,
    diagX = diagX
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
#' @param iter_selection Number of post-burn-in samples for model selection.
#' @param candidates Candidate models and their initial a,b values.
#' @param woodbury Logical, use Woodbury identity in beta updates.
#' @param delta1,delta2,delta3 Convergence tolerances.
#' @param window_size Size of the convergence window.
#' @param diagX Logical, assume diagonal X.
#' @return A list containing winning_params and winner_name.
#' @export
initialize_adaptive <- function(X, y,
                                omega_init_guess,
                                sigmaSq_init_guess,
                                iter_pre_opt = 100,
                                iter_selection = 5000,
                                candidates = list(
                                  horseshoe = list(a = 0.5, b = 0.5),
                                  normal_gamma = list(a = 0.5, b = 5.0),
                                  studentt = list(a = 5.0, b = 0.5)
                                ),
                                woodbury = TRUE,
                                delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                                window_size = 5, 
                                diagX = FALSE) {
  cat("Stage 1: Pre-optimizing each candidate model via Laplace-ECM...\n")

  pre_optimized_params <- list()
  pre_optimized_states <- list()
  engine_args <- list(
    X = X, y = y,
    max_iter = iter_pre_opt,
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
        update_a = FALSE,
        update_b = FALSE,
        update_omega = TRUE,
        update_sigmaSq = TRUE
      ),
      engine_args
    ))
    pre_optimized_params[[name]] <- opt_res$params
    pre_optimized_states[[name]] <- opt_res$final_state
  }

  cat("Stage 2: Fast pre-sampling via Laplace Approximation ...\n")
  presampled_betas <- list()
  n <- nrow(X)

  for (name in names(candidates)) {
    starting <- pre_optimized_states[[name]]
    params <- pre_optimized_params[[name]]
    p <- length(starting$beta0)
    if (is.null(starting$nu0) || is.null(starting$lambda0)) {
      stop("Stage 2 requires final nu0 and lambda0 from eb_run_em_engine().")
    }
    sigmaSq_current <- eb_positive_finite(params$sigmaSq, 1e-16)
    omega_current <- eb_positive_finite(params$omega, 1e-16)
    raw_d <- omega_current * starting$nu0 / starting$lambda0
    d_current <- pmax(ifelse(is.finite(raw_d), raw_d, .Machine$double.xmin),
                      .Machine$double.xmin)
    
    if (woodbury && !diagX) {
      XD <- sweep(X, 2, d_current, FUN = "*")
      M <- XD %*% t(X)
      diag(M) <- diag(M) + sigmaSq_current
      M_inv <- eb_chol2inv_pd(M)
      U <- matrix(rnorm(p * iter_selection, mean = 0, sd = sqrt(d_current)), nrow = p, ncol = iter_selection)
      Delta <- matrix(rnorm(n * iter_selection, mean = 0, sd = 1), nrow = n, ncol = iter_selection)
      Y_mat <- matrix(y, nrow = n, ncol = iter_selection)
      RHS <- Y_mat - (X %*% U) - sqrt(sigmaSq_current) * Delta
      W <- M_inv %*% RHS
      beta_draws <- t(U + t(XD) %*% W)
    } else {
      if (diagX) {
        beta_draws <- matrix(
          rnorm(iter_selection * p, 
                mean = rep(starting$beta0, each = iter_selection), 
                sd = rep(sqrt(pmax(starting$v_beta, 0)), each = iter_selection)),
          nrow = iter_selection, ncol = p
        )
      } else {
        XtX <- crossprod(X)
        precision_beta <- XtX / sigmaSq_current + diag(1 / d_current, p)
        sigma_beta <- eb_chol2inv_pd(precision_beta)
        beta_draws <- MASS::mvrnorm(n = iter_selection, mu = starting$beta0, Sigma = sigma_beta)
      }
    }
    
    presampled_betas[[name]] <- beta_draws
  }

  cat("Stage 3: Starting iterative model selection ...\n")
  model_names <- names(candidates)
  prob_matrix <- matrix(NA, nrow = iter_selection, ncol = length(candidates))
  current_model_name <- sample(model_names, 1)
  current_beta <- presampled_betas[[current_model_name]][1, ]

  for (j in seq_len(iter_selection)) {
    log_weights <- sapply(model_names, function(name) {
      eb_calculate_marginal_loglik_beta(
        beta_vec = current_beta,
        model_params = pre_optimized_params[[name]]
      )
    })

    max_log <- max(log_weights, na.rm = TRUE)
    if (!is.finite(max_log)) {
      probs <- rep(1 / length(model_names), length(model_names))
    } else {
      weights <- exp(log_weights - max_log)
      probs <- weights / sum(weights, na.rm = TRUE)
    }
    
    current_model_name <- sample(model_names, 1, prob = probs)
    sample_row_idx <- sample(1:iter_selection, 1)
    
    current_beta <- presampled_betas[[current_model_name]][sample_row_idx, ]

    prob_matrix[j, ] <- probs
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



eb_initial_values <- function(X, y, option = "ridge", epsilon = 1e-6) {
  n <- nrow(X)

  if (option == "ridge") {
    cv_ridge <- glmnet::cv.glmnet(X, y, alpha = 0, intercept = FALSE)
    y_hat <- as.vector(predict(cv_ridge, s = "lambda.min", newx = X))
    sigmaSq_hat <- max(epsilon, var(as.vector(y) - y_hat, na.rm = TRUE))
    
  } else if (option == "olasso") {
    ol_fit <- natural::olasso(X, y, intercept = FALSE)
    sigmaSq_hat <- max(epsilon, ol_fit$sig_obj_1)
    
  } else {
    stop("Invalid option. Please choose 'ridge' or 'olasso'.")
  }

  y_sq_sum <- sum(y^2)
  trace_XX <- sum(X^2)
  # From: 1/n * y'y = 1/n * trace(X'X) * omega + sigmaSq
  # omega = (y'y - n * sigmaSq) / trace(X'X)
  numerator <- y_sq_sum - (n * sigmaSq_hat)
  if (trace_XX <= 0 || numerator <= 0) {
    omega_hat <- epsilon
  } else {
    omega_hat <- max(epsilon, numerator / trace_XX)
  }
  return(list(sigmaSq = sigmaSq_hat, omega = omega_hat))
}


eb_run_em_engine <- function(X, y,
                             a_init, b_init, omega_init, sigmaSq_init,
                             max_iter,
                             delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3, window_size = 5,
                             woodbury = TRUE,
                             update_a = FALSE,
                             update_b = FALSE,
                             update_omega = TRUE,
                             update_sigmaSq = TRUE,
                             diagX = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  y <- as.vector(y)
  a_vec <- c(a_init)
  b_vec <- c(b_init)
  omega_vec <- c(omega_init)
  sigmaSq_vec <- c(sigmaSq_init)
  beta_diff_vec <- numeric(0)
  beta0 <- rep(0, p)
  nu0 <- rep(1, p)
  lambda0 <- rep(1, p)
  XtX <- if (diagX) NULL else crossprod(X)
  Xty <- if (diagX) NULL else crossprod(X, y)

  for (k in 1:max_iter) {
    a_current <- a_vec[k]
    b_current <- b_vec[k]
    omega_current <- eb_positive_finite(omega_vec[k], 1e-16)
    sigmaSq_current <- eb_positive_finite(sigmaSq_vec[k], 1e-16)

    raw_d <- omega_current * nu0 / lambda0
    d_current <- pmax(ifelse(is.finite(raw_d), raw_d, .Machine$double.xmin),
                      .Machine$double.xmin)
    if (!diagX) {
      if (woodbury) {
        XD <- sweep(X, 2, d_current, FUN = "*") 
        M <- XD %*% t(X)                       
        diag(M) <- diag(M) + sigmaSq_current    
        M_inv <- eb_chol2inv_pd(M)

        w <- M_inv %*% y
        beta_mode <- as.vector(t(XD) %*% w)
        MinvX <- M_inv %*% X
        diag_XtMinvX <- colSums(X * MinvX) 
        v_beta <- pmax(d_current - (d_current^2) * diag_XtMinvX, 0)

        trace_term <- sigmaSq_current * (n - sigmaSq_current * sum(diag(M_inv)))
        residuals_mode <- as.vector(y - X %*% beta_mode)
      } else {
        precision_beta <- XtX / sigmaSq_current + diag(1 / d_current, p)
        sigma_beta <- eb_chol2inv_pd(precision_beta)
        beta_mode <- as.vector(sigma_beta %*% (Xty / sigmaSq_current))
        v_beta <- pmax(diag(sigma_beta), 0)
        residuals_mode <- y - as.vector(X %*% beta_mode)
        trace_term <- sum(XtX * sigma_beta)
      }
    } else {
      x_diag <- diag(X)
      precision_beta_diag <- x_diag^2 / sigmaSq_current + 1 / d_current
      v_beta <- pmax(1 / precision_beta_diag, 0)
      beta_mode <- v_beta * x_diag * y / sigmaSq_current
      residuals_mode <- y - x_diag * beta_mode
      trace_term <- sum(x_diag^2 * v_beta)
    }

    p_gig <- a_current - 0.5
    a_gig <- 2 * a_current
    b_gig <- pmax(beta_mode^2 * lambda0 / omega_current, 0)
    gig_root <- sqrt(p_gig^2 + a_gig * b_gig)
    nu_mode <- pmax((p_gig + gig_root) / a_gig, .Machine$double.xmin)

    lambda_rate <- b_current + beta_mode^2 / (2 * omega_current * nu_mode)
    lambda_mode <- pmax((b_current + 0.5) / lambda_rate, .Machine$double.xmin)

    a_new <- if (update_a) {
      max(1e-6, eb_idigamma(mean(log(nu_mode)) + log(a_current)))
    } else {
      a_current
    }
    b_new <- if (update_b) {
      max(1e-6, eb_idigamma(mean(log(lambda_mode)) + log(b_current)))
    } else {
      b_current
    }
    sigmaSq_new <- if (update_sigmaSq) {
      max((sum(residuals_mode^2) + trace_term) / n, 1e-16)
    } else {
      sigmaSq_current
    }
    omega_new <- if (update_omega) {
      eb_laplace_ecm_omega(
        beta_mode = beta_mode,
        v_beta = v_beta,
        d_current = d_current,
        omega_current = omega_current
      )
    } else {
      omega_current
    }

    beta_diff_vec <- c(beta_diff_vec, mean((beta_mode - beta0)^2))
    beta0 <- beta_mode
    nu0 <- nu_mode
    lambda0 <- lambda_mode
    a_vec <- c(a_vec, a_new)
    b_vec <- c(b_vec, b_new)
    sigmaSq_vec <- c(sigmaSq_vec, sigmaSq_new)
    omega_vec <- c(omega_vec, omega_new)

    converged <- eb_isConverged_hyper(a_vec, b_vec, sigmaSq_vec, omega_vec, delta1 = delta1, delta2 = delta2, window_size = window_size)

    if (converged && length(beta_diff_vec) >= window_size &&
        max(tail(beta_diff_vec, window_size), na.rm = TRUE) < delta3) {
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
    final_state = list(beta0 = beta0, v_beta = v_beta,
                       nu0 = nu0, lambda0 = lambda0),
    trajectories = list(a_traj = a_vec, b_traj = b_vec,
                        sigmaSq_traj = sigmaSq_vec,
                        omega_traj = omega_vec,
                        beta_diff_traj = beta_diff_vec)
  )
}



eb_laplace_ecm_omega <- function(beta_mode, v_beta, d_current, omega_current) {
  # Numerically stabilized conditional expectation.
  ez2 <- omega_current * (beta_mode^2 + v_beta) / d_current
  ez2 <- ez2[is.finite(ez2)]
  if (length(ez2) == 0) {
    return(omega_current)
  }
  return(max(mean(ez2), 1e-16))
}

eb_positive_finite <- function(x, lower) {
  if (!is.finite(x) || x <= lower) {
    return(lower)
  }
  x
}

eb_chol2inv_pd <- function(mat, jitter = 1e-10, max_tries = 6) {
  mat <- (mat + t(mat)) / 2
  for (i in 0:max_tries) {
    bump <- if (i == 0) 0 else jitter * 10^(i - 1)
    candidate <- mat
    diag(candidate) <- diag(candidate) + bump
    chol_candidate <- tryCatch(chol(candidate), error = function(e) NULL)
    if (!is.null(chol_candidate)) {
      return(chol2inv(chol_candidate))
    }
  }
  stop("Matrix is not positive definite after jitter stabilization.")
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
  const <- gamma(0.5 + b) * gamma(a + b) / (gamma(a) * gamma(b) * sqrt(2 * pi * phi))
  U_a <- 0.5 + b
  U_b <- 1.5 - a
  
  vapply(beta_vec, function(beta_i) {
    z_val <- max(beta_i^2 / (2 * phi), .Machine$double.xmin)
    
    U_val <- tryCatch({
      gsl::hyperg_U(U_a, U_b, z_val)
    }, error = function(e) NA_real_)
    
    if (is.na(U_val)) {
      return(NA_real_)
    }
    pmin(const * U_val, .Machine$double.xmax)
  }, numeric(1))
}
