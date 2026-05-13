#' @import mvtnorm
#' @import glmnet
#' @param mainEM_burnin number of burn-in samples in each iteration in the main EM.
#' @param mainEM_samples number of samples used for estimation in each iteration in the main EM
#' @param output_burnin number of burn-in samples needed for the draw after convergence
#' @param output_samples number of samples produced after convergence
#' @param delta1 numeric, a small number added to avoid division by zero
#' @param delta2 numeric, error bound for convergence in isConverged_hyper
#' @param delta3 numeric, error bound for convergence in isConverged_rho
#' @param window_size numeric, size of the sliding window
#' @param converge_rule "rho" or "hyper", specifying what to monitor to decide the convergence
#' @param mainEM_iter numeric, max number of iterations in the main EM.
#' @param IS_period number of iterations between each resampling step in importance sampling. Set to 1 by default to disable importance sampling
#' @param option "augmented" or "naive", specifying how we would like to sample nu
#' @param a.start starting point of a for EM-within-Gibbs
#' @param a a,b,omega,sigmaSq, set to some specific values to if we would like to fix them, otherwise set to NULL
#' @param adapt_init TRUE if we want to use the initialize_adaptive function to find starting points
#' @return list of posterior samples
#' @export
EM_GibbsTPB_adapt = function(X, y,
                             iter_pre_opt, pre_opt_burnin, pre_opt_samples,
                             iter_burnin_selection, iter_selection,
                             mainEM_burnin, mainEM_samples,
                             output_burnin, output_samples,
                             mainEM_iter = 500,
                             IS_period = 1, option = "augmented",
                             woodbury = FALSE, approx = FALSE,
                             delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                             window_size = 5, converge_rule = "hyper",
                             a = NULL, b = NULL, omega = NULL, sigmaSq = NULL,
                             a.start = 0.5, b.start = 0.5, omega.start = 1, sigmaSq.start = 1,
                             adapt_init = TRUE, diagX = FALSE){
  
  n = nrow(X)
  p = ncol(X)
  
  
  # --- 1. INITIALIZATION STEP ---
  start_params <- list(a = a.start, b = b.start, omega = omega.start, sigmaSq = sigmaSq.start)
  
  if (adapt_init) {
    if (!is.null(a) || !is.null(b) || !is.null(omega) || !is.null(sigmaSq)) {
      cat("Hyperparameters are fixed by user, skipping adaptive initialization.\n")
    } else {
      start_model <- initialize_adaptive(
        X = X, y = y,
        omega_init_guess = omega.start, sigmaSq_init_guess = sigmaSq.start,
        iter_pre_opt = iter_pre_opt, pre_opt_burnin = pre_opt_burnin, pre_opt_samples = pre_opt_samples,
        iter_burnin_selection = iter_burnin_selection, iter_selection = iter_selection,
        woodbury = woodbury, IS_period = IS_period,
        delta1 = delta1, delta2 = delta2, delta3 = delta3,
        window_size = window_size, converge_rule = converge_rule, 
        approx = approx, option = option, diagX = diagX
      )
      start_params = start_model$winning_params
    }
  }
  
  # --- 2. MAIN EM RUN STEP ---
  cat("\n--- Starting Main EM Run ---\n")
  final_em_results <- run_em_engine(
    X = X, y = y,
    a_init = start_params$a, b_init = start_params$b, 
    omega_init = start_params$omega, sigmaSq_init = start_params$sigmaSq,
    null.a = is.null(a), null.b = is.null(b), 
    null.omega = is.null(omega), null.sigmaSq = is.null(sigmaSq),
    max_iter = mainEM_iter, iter_burnin = mainEM_burnin, iter_samples = mainEM_samples, 
    delta1 = delta1, delta2 = delta2, delta3 = delta3, 
    window_size = window_size, converge_rule = converge_rule,
    woodbury = woodbury, approx = approx,
    IS_period = IS_period, option = option, diagX = diagX
  )
  
  # --- 3. FINAL POSTERIOR SAMPLING STEP ---
  cat("\n--- Generating Final Posterior Samples ---\n")
  final_params <- final_em_results$params
  final_sampler_state <- final_em_results$final_state
  
  output <- getsamples.emp(
    num = output_samples, X = X, y = y,
    a = final_params$a, b = final_params$b, omega = final_params$omega, sigmaSq = final_params$sigmaSq,
    beta0 = final_sampler_state$beta0, nu0 = final_sampler_state$nu0, lambda0 = final_sampler_state$lambda0,
    burn = output_burnin, xi0 = final_sampler_state$xi0, woodbury = woodbury, approx = approx, diagX = diagX
  )
  
  res <- list(
    trajectories = final_em_results$trajectories,
    final_params = final_params,
    beta = output$beta
  )
  if (adapt_init) {
    res$start_model = start_model$winner_name
  }
  
  return(res)
}