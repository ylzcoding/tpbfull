#' Performs the Data-Adaptive Initialization Strategy via Candidate Model Competition.
#'
#' This function automates the selection of starting points for the main EM_GibbsTPB
#' algorithm by running a competition between several candidate models (special cases
#' of the TPB prior). It selects the candidate that best fits the data and returns
#' its parameters as the starting points.
#' 
#' @param X, y Model matrix and response vector.
#' @param omega_init_guess, sigmaSq_init_guess Initial guesses for omega and sigmaSq, typically derived from a fast Ridge fit.
#' @param iter_pre_opt Iterations for the pre-optimization EM stage for each candidate.
#' @param pre_opt_burnin, pre_opt_samples burn-in and sample size for the Gibbs sampler within each pre-optimization EM step.
#' @param iter_burnin_selection Number of burn-in iterations for the main model selection MCMC chain.
#' @param iter_selection Number of post-burn-in samples for model selection.
#' @param candidates A list defining the candidate models and their fixed (a, b) values.
#' @param ... Other arguments to be passed to run_em_engine and getsamples.emp.
#' @return A list containing the winning parameters and the name of the winning model.
#' @export

initialize_adaptive <- function(X, y,
                                omega_init_guess,
                                sigmaSq_init_guess,
                                iter_pre_opt = 200,
                                pre_opt_burnin = 1000, pre_opt_samples = 1000,
                                iter_burnin_selection = 0,
                                iter_selection = 5000,
                                candidates = list(
                                  horseshoe = list(a = 0.5, b = 0.5),
                                  sb = list(a = 1.0, b = 0.5),
                                  normal_gamma = list(a = 0.5, b = 20.0),
                                  studentt = list(a = 20.0, b = 1.0)
                                ),
                                woodbury = TRUE, IS_period = 1,
                                delta1 = 1e-6, delta2 = 1e-3, delta3 = 1e-3,
                                window_size = 5, converge_rule = "hyper",
                                approx = FALSE, option = "augmented", diagX = FALSE) {
  
  
  cat("Stage 1: Pre-optimizing each candidate model...\n")
  
  # iter_pre_opt, pre_opt_burnin, and pre_opt_samples are parameters for pre-optimizing the candidates. 
  # For each candidate, we fix a and b, and use omega_init_guess and sigmaSq_init_guess as starting values 
  # to run the pre-optimization using em-within-gibbs, so that omega and sigmaSq converge to some better values.
  
  # iter_pre_opt: Iterations for the pre-optimization EM stage for each candidate
  # pre_opt_burnin, pre_opt_samples: burn-in and post burnin sample size for the Gibbs sampler within each pre-optimization EM step.
  
  pre_optimized_params <- list()
  pre_optimized_states <- list()
  
  engine_args <- list(
    X = X, y = y, 
    max_iter = iter_pre_opt, 
    iter_burnin = pre_opt_burnin,
    iter_samples = pre_opt_samples,
    delta1 = delta1, delta2 = delta2, delta3 = delta3,
    window_size = window_size, converge_rule = converge_rule,
    woodbury = woodbury, approx = approx,
    IS_period = IS_period, option = option, diagX = diagX)
  
  for (name in names(candidates)) {
    cat("Optimizing candidate model:", name, "...\n")
    candidate_ab <- candidates[[name]]
    specific_args <- list(
      a_init = candidate_ab$a, b_init = candidate_ab$b,
      omega_init = omega_init_guess, sigmaSq_init = sigmaSq_init_guess,
      null.a = FALSE, null.b = FALSE,
      null.omega = TRUE, null.sigmaSq = TRUE
    )
    opt_res <- do.call(run_em_engine, c(specific_args, engine_args))
    pre_optimized_params[[name]] <- opt_res$params
    pre_optimized_states[[name]] <- opt_res$final_state
    
    # For each candidate, at the end of pre-optimization we store its pre_optimized_params, 
    # along with the beta, lambda, and nu samples obtained in the final iteration of the EM algorithm.
  }
  
  # get pre-sampled beta, lambda, nu, and xi (#pre-sampled beta = iter_selection)
  
  cat("Stage 2: pre-sampling ...\n")
  
  presampled_betas <- list()
  presampled_nus <- list()
  presampled_lambdas <- list()
  
  for (name in names(candidates)) {
    params <- pre_optimized_params[[name]]
    starting <- pre_optimized_states[[name]]
    presamples <- getsamples.emp(num = iter_selection, X = X, y = y,
                                 a = params$a, b = params$b, 
                                 omega = params$omega, sigmaSq = params$sigmaSq,
                                 beta0 = starting$beta0, nu0 = starting$nu0, 
                                 lambda0 = starting$lambda0, burn = 0,
                                 xi0 = starting$xi0, woodbury=woodbury, approx=approx, diagX=diagX)
    presampled_betas[[name]] <- presamples$beta
    presampled_nus[[name]] <- presamples$nu
    presampled_lambdas[[name]] <- presamples$lambda
  }
  
  
  cat("Stage 3: Starting iterative model selection ...\n")
  
  # iter_burnin_selection: number of burn-in iterations for the main model selection MCMC chain (if applicable)
  # iter_selection: number of post-burnin samples for model selection.
  
  
  num_candidates <- length(candidates)
  model_names <- names(candidates)
  
  
  # initialize a probability matrix of size iter_selection × num_candidates
  # to store the multinomial probabilities and compute the average at the end
  prob_matrix <- matrix(NA, nrow = iter_selection, ncol = num_candidates)
  # randomly generate a delta at the very beginning
  current_model_name <- sample(model_names, 1)
  
  current_beta <- presampled_betas[[current_model_name]][1, ]
  current_nu <- presampled_nus[[current_model_name]][1, ]
  current_lambda <- presampled_lambdas[[current_model_name]][1, ]
  
  
  # the total number of iterations
  iter_total <- iter_burnin_selection + iter_selection
  
  for (j in 1:iter_total) {
    # Given current beta, lambda and nu, we compute the multinomial probabilities
    # This step is to update model indicator variable, delta
    log_weights <- sapply(model_names, function(name) {
      params <- pre_optimized_params[[name]]
      calculate_marginal_loglik_beta(beta_vec = current_beta, model_params = params)
    })
    
    max_log_weight <- max(log_weights)
    weights <- exp(log_weights - max_log_weight)
    probs <- weights / sum(weights, na.rm = TRUE)
    
    # resample delta from the resulting multinomial distribution
    current_model_name <- sample(model_names, 1, prob = probs)
    
    # next, based on the current candidate model, we retrieve a new beta, nu, and lambda from the presampled array 
    sample_row_idx <- sample(1:iter_selection, 1)
    current_beta <- presampled_betas[[current_model_name]][sample_row_idx, ]
    current_nu <- presampled_nus[[current_model_name]][sample_row_idx, ]
    current_lambda <- presampled_lambdas[[current_model_name]][sample_row_idx, ]
    
    # Only store choices after the burn-in period
    if (j > iter_burnin_selection) {
      storage_idx <- j - iter_burnin_selection
      prob_matrix[storage_idx, ] <- probs
    }
  }
  
  avg_probs <- colMeans(prob_matrix)
  names(avg_probs) <- model_names
  winner_name <- names(which.max(avg_probs))
  winning_params <- pre_optimized_params[[winner_name]]
  
  
  return(
    list(winning_params = winning_params,
         winner_name = winner_name)
  )
}
  
  
  
  
  
  
  
                                  
