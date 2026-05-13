#' The core EM-within-Gibbs engine.
#' This function runs the iterative EM algorithm until convergence or max_iter is reached.
#' It can handle fixed or updating hyperparameters based on the 'null.*' arguments.
#'
#' @param iter_burnin number of burn-in samples in each iteration
#' @param iter_samples number of samples used for estimation in each iteration
#' @param a_init, b_init, omega_init, sigmaSq_init Initial values for hyperparameters.
#' @param null.a, null.b, null.omega, null.sigmaSq Booleans indicating which parameters to update (TRUE) vs. keep fixed (FALSE).
#' @param max_iter Maximum number of iterations for this run.
#' @param IS_period The frequency of resampling vs. using IS
#' @param option "augmented" or "naive", specifying how we would like to sample nu
#' @param converge_rule "rho" or "hyper", specifying what to monitor to decide the convergence
#' @param window_size numeric, size of the sliding window
#' @return A list containing the converged parameters (params) and the final state of the samplers (final_state).
run_em_engine <- function(X, y,
                          a_init, b_init, omega_init, sigmaSq_init,
                          null.a, null.b, null.omega, null.sigmaSq,
                          max_iter, iter_burnin, iter_samples,
                          delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                          window_size = 5, converge_rule = "hyper",
                          woodbury = FALSE, approx = FALSE,
                          IS_period = 1, option = "augmented", diagX = FALSE
                          ) {
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  # Initialize vectors to record trajectories
  a_vec <- c(a_init)
  b_vec <- c(b_init)
  omega_vec <- c(omega_init)
  sigmaSq_vec <- c(sigmaSq_init)

  beta_diff_vec <- numeric(0) 
  last_beta_hat <- NULL 
  
  # Initialize latent variables
  beta0 <- rnorm(p, mean = 0, sd = sqrt(sigmaSq_init))
  nu0 <- rgamma(p, 1, 1)
  lambda0 <- rgamma(p, 1, 1)
  if (option == "augmented") {
    xi0 = a_init * sqrt(exp(log(beta0^2) + log(lambda0) - log(omega_init))/2)
  } else {
    xi0 = NULL
  }

  resample_iter <- 0
  
  # --- Start of the EM Loop ---
  for (k in 1:max_iter) {
    
    current_beta_hat <- NULL
    
    if ((k - 1) %% IS_period == 0) {
      resample_iter <- k
      samples <- getsamples.emp(
        num = iter_samples, X = X, y = y,
        a = a_vec[k], b = b_vec[k], omega = omega_vec[k], sigmaSq = sigmaSq_vec[k],
        beta0 = beta0, nu0 = nu0, lambda0 = lambda0, burn = iter_burnin, xi0 = xi0,
        woodbury = woodbury, approx = approx, diagX = diagX
      )
      
      # Regular M-step, update parameters based on the 'null.*' flags
      a_new <- if (null.a) M.step(sample_matrix = samples$nu, param_current = a_vec[k]) else a_vec[k]
      b_new <- if (null.b) M.step(sample_matrix = samples$lambda, param_current = b_vec[k]) else b_vec[k]
      sigmaSq_new <- if (null.sigmaSq) M.step_sigmaSq(samples$beta, X, y, nrow(X), diagX) else sigmaSq_vec[k]
      omega_new <- if (null.omega) M.step_phi(samples$beta, samples$lambda, samples$nu) else omega_vec[k]
      
      # Update initial values for the next full E-step
      beta0 <- samples$beta[iter_samples, ]
      nu0 <- samples$nu[iter_samples, ]
      lambda0 <- samples$lambda[iter_samples, ]
      xi0 <- samples$xi[iter_samples, ]

      current_beta_hat <- colMeans(samples$beta)
      
    } else {
      # M-step with Importance Sampling
      weights_vec <- IS_weights(y = y, X = X,
                                beta_mat = samples$beta, lambda_mat = samples$lambda, nu_mat = samples$nu,
                                a_cur = a_vec[k], b_cur = b_vec[k], omega_cur = omega_vec[k], sigmaSq_cur = sigmaSq_vec[k],
                                a0 = a_vec[resample_iter], b0 = b_vec[resample_iter], omega0 = omega_vec[resample_iter], sigmaSq0 = sigmaSq_vec[resample_iter])
      
      
      a_new <- if (null.a) M.step.IS(weights_vec, samples$nu * a_vec[resample_iter]) else a_vec[k]
      b_new <- if (null.b) M.step.IS(weights_vec, samples$lambda * b_vec[resample_iter]) else b_vec[k]
      sigmaSq_new <- if (null.sigmaSq) M.step_sigmaSq.IS(weights_vec, samples$beta, X, y, n, diagX) else sigmaSq_vec[k]
      omega_new <- if (null.omega) M.step_omega.IS(weights_vec, samples$beta, samples$lambda, samples$nu) else omega_vec[k]

      sum_weights <- sum(weights_vec)
      if (sum_weights > 0) {
        current_beta_hat <- colSums(samples$beta * weights_vec) / sum_weights
      } else {
        current_beta_hat <- colMeans(samples$beta) # fallback
      }
    }
      
    if (!is.null(last_beta_hat)) {
      diff_mse <- mean((current_beta_hat - last_beta_hat)^2)
      beta_diff_vec <- c(beta_diff_vec, diff_mse)
    } else {
      beta_diff_vec <- c(beta_diff_vec, NA)
    }
    
    # update last_beta_hat for next iteration
    last_beta_hat <- current_beta_hat

    a_vec <- c(a_vec, a_new)
    b_vec <- c(b_vec, b_new)
    sigmaSq_vec <- c(sigmaSq_vec, sigmaSq_new)
    omega_vec <- c(omega_vec, omega_new)
    
    # Convergence Check
    if (converge_rule == "hyper") {
      if (isConverged_hyper(a_vec, b_vec, sigmaSq_vec, omega_vec, delta1, delta2, window_size)) {
            cat("Engine converged after", k, "iterations.\n")
            break
        }
      } else {
        if (isConverged_rho(a_vec, b_vec, omega_vec, delta3, window_size)) {
          cat("Engine converged after", k, "iterations.\n")
          break
            }
      }
    cat(k, "-th iteration completed.\n")
    cat("a=", a_vec[k+1], "\n")
    cat("b=", b_vec[k+1], "\n")
    cat("omega=", omega_vec[k+1], "\n")
    cat("sigmaSq=", sigmaSq_vec[k+1], "\n")
  }
  
  # --- End of the EM Loop ---
  
  l <- length(a_vec)
  final_params <- list(
    a = a_vec[l], b = b_vec[l], omega = omega_vec[l], sigmaSq = sigmaSq_vec[l]
  )
  final_state <- list(beta0 = beta0, nu0 = nu0, lambda0 = lambda0, xi0 = xi0)
  return(
    list(params = final_params, 
         final_state = final_state,
         trajectories = list(a_traj = a_vec, b_traj = b_vec, sigmaSq_traj = sigmaSq_vec, omega_traj = omega_vec, beta_diff_traj = beta_diff_vec)
    )
  )
}
