#' TPB Full Pipeline with Optional EB Prior Elicitation and Multi-Chain MCMC
#'
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param num_chains Integer, number of independent MCMC chains to run
#' @param auto_find Logical, run empirical Bayes model competition once before MCMC
#' @param hyper_params List of hyperparameters for the fully Bayesian sampler
#' @param num_iter Integer, total number of MCMC iterations across all chains
#' @param num_warmup Integer, total number of warm-up iterations across all chains
#' @param thinning Integer, thinning interval for saved posterior samples
#' @param iter_pre_opt Iterations for the model competition pre-optimization EM stage for each candidate
#' @param pre_opt_burnin Burn-in for each model competition pre-optimization EM step
#' @param pre_opt_samples Sample size for each model competition pre-optimization EM step
#' @param iter_selection Number of post-burn-in samples for model selection
#' @param woodbury Logical, use Woodbury identity in both model competition and full Gibbs sampling
#' @param ... Additional arguments passed to fullGibbs()
#' @return A list containing combined samples, per-chain outputs, diagnostics,
#'   hyperparameters used, and optional EB competition metadata.
#' @export
tpb_full_pipeline <- function(X, y,
                              num_chains = 4,
                              auto_find = TRUE,
                              hyper_params = list(prior_type_a = "gamma",
                                                  prior_type_b = "gamma",
                                                  s_a = 1.5, r_a = 1,
                                                  s_b = 1.5, r_b = 1,
                                                  scale_a = 1,
                                                  scale_b = 1,
                                                  scale_phi = 1),
                              num_iter = 100000,
                              num_warmup = 25000,
                              thinning = 1,
                              iter_pre_opt = 100,
                              pre_opt_burnin = 1000,
                              pre_opt_samples = 1000,
                              iter_selection = 5000,
                              woodbury = TRUE,
                              ...) {


  num_posterior <- num_iter - num_warmup
  warmup_per_chain <- num_warmup %/% num_chains
  posterior_per_chain <- num_posterior %/% num_chains
  iter_per_chain <- warmup_per_chain + posterior_per_chain

  if (posterior_per_chain < 1L) {
    stop("num_iter - num_warmup must be at least num_chains so each chain saves at least one posterior draw.")
  }

  final_hyper_params <- modifyList(
    list(prior_type_a = "gamma", prior_type_b = "gamma",
         s_a = 1.5, r_a = 1, s_b = 1.5, r_b = 1,
         scale_a = 1, scale_b = 1, scale_phi = 1),
    hyper_params
  )

  competition_result <- NULL
  winning_modes <- NULL
  if (isTRUE(auto_find)) {
    competition_result <- run_model_competition(
      X = X,
      y = y,
      iter_pre_opt = iter_pre_opt,
      pre_opt_burnin = pre_opt_burnin,
      pre_opt_samples = pre_opt_samples,
      iter_selection = iter_selection,
      woodbury = woodbury
    )
    winning_modes <- list(
      a = competition_result$winner$a,
      b = competition_result$winner$b
    )
    final_hyper_params <- modifyList(
      final_hyper_params,
      list(
        prior_type_a = "gamma",
        prior_type_b = "gamma",
        s_a = winning_modes$a + 1,
        r_a = 1,
        s_b = winning_modes$b + 1,
        r_b = 1
      )
    )
  }

  chain_results <- vector("list", num_chains)
  chains_beta <- vector("list", num_chains)
  chains_scalars <- vector("list", num_chains)
  chains_loglik <- vector("list", num_chains)
  chains_logpost <- vector("list", num_chains)
  chains_accept <- vector("list", num_chains)
  chains_cov <- vector("list", num_chains)
  chains_scale <- vector("list", num_chains)

  for (i in seq_len(num_chains)) {
    cat(sprintf("Running fullGibbs chain %d / %d ...\n", i, num_chains))
    chain_res <- fullGibbs(
      X = X,
      y = y,
      num_output = posterior_per_chain * thinning,
      num_burnin = warmup_per_chain * thinning,
      thin = thinning,
      woodbury = woodbury,
      hyper_params = final_hyper_params,
      ...
    )

    chain_results[[i]] <- chain_res
    chains_beta[[i]] <- chain_res$samples$beta
    chains_scalars[[i]] <- chain_res$samples$scalars
    chains_loglik[[i]] <- chain_res$samples$beta_loglik
    chains_logpost[[i]] <- chain_res$samples$total_logpost
    chains_accept[[i]] <- unlist(chain_res$acceptance_rates[c("a", "b", "phi")])
    chains_cov[[i]] <- chain_res$final_proposal_cov
    chains_scale[[i]] <- chain_res$final_scale_factor
  }

  combined_beta <- do.call(rbind, chains_beta)
  combined_scalars <- do.call(rbind, chains_scalars)
  combined_loglik <- unlist(chains_loglik, use.names = FALSE)
  combined_logpost <- unlist(chains_logpost, use.names = FALSE)

  accept_matrix <- do.call(rbind, chains_accept)
  avg_acceptance <- colMeans(accept_matrix)
  diagnostics <- tpb_mcmc_diagnostics(chains_scalars)

  has_adaptive_tuning <- any(vapply(chains_cov, Negate(is.null), logical(1)))
  avg_cov <- if (has_adaptive_tuning) {
    Reduce("+", chains_cov[!vapply(chains_cov, is.null, logical(1))]) /
      sum(!vapply(chains_cov, is.null, logical(1)))
  } else {
    NULL
  }
  avg_scale <- if (any(vapply(chains_scale, Negate(is.null), logical(1)))) {
    mean(unlist(chains_scale[!vapply(chains_scale, is.null, logical(1))]))
  } else {
    NULL
  }

  list(
    samples = list(
      beta = combined_beta,
      scalars = combined_scalars,
      beta_loglik = combined_loglik,
      total_logpost = combined_logpost
    ),
    chains = list(
      results = chain_results,
      beta = chains_beta,
      scalars = chains_scalars,
      beta_loglik = chains_loglik,
      total_logpost = chains_logpost,
      acceptance_rates = chains_accept,
      final_proposal_cov = chains_cov,
      final_scale_factor = chains_scale
    ),
    scalar_traces = list(
      a_samples = as.vector(combined_scalars[, "a"]),
      b_samples = as.vector(combined_scalars[, "b"]),
      phi_samples = as.vector(combined_scalars[, "phi"])
    ),
    acceptance_rates = list(
      per_chain = accept_matrix,
      mean = avg_acceptance
    ),
    adaptive_tuning = list(
      per_chain_cov = chains_cov,
      per_chain_scale = chains_scale,
      mean_cov = avg_cov,
      mean_scale = avg_scale
    ),
    diagnostics = diagnostics,
    mcmc_schedule = list(
      num_chains = num_chains,
      num_iter_total = num_iter,
      num_warmup_total = num_warmup,
      num_posterior_total = num_posterior,
      thinning = thinning,
      iter_per_chain = iter_per_chain,
      warmup_per_chain = warmup_per_chain,
      posterior_per_chain = posterior_per_chain
    ),
    hyper_params = final_hyper_params,
    auto_find = list(
      enabled = isTRUE(auto_find),
      competition_result = competition_result,
      winning_modes = winning_modes
    )
  )
}

tpb_mcmc_diagnostics <- function(chains_scalars) {
  if (length(chains_scalars) == 0) {
    return(list(
      rhats = numeric(0),
      ess = numeric(0),
      max_rhat = NA_real_,
      min_ess = NA_real_
    ))
  }

  scalar_names <- colnames(chains_scalars[[1]])
  ess_vals <- tryCatch({
    if (length(chains_scalars) == 1) {
      coda::effectiveSize(coda::mcmc(chains_scalars[[1]]))
    } else {
      coda::effectiveSize(coda::mcmc.list(lapply(chains_scalars, coda::mcmc)))
    }
  }, error = function(e) {
    stats::setNames(rep(NA_real_, length(scalar_names)), scalar_names)
  })

  rhats <- if (length(chains_scalars) < 2) {
    stats::setNames(rep(NA_real_, length(scalar_names)), scalar_names)
  } else {
    tryCatch({
      coda::gelman.diag(
        coda::mcmc.list(lapply(chains_scalars, coda::mcmc)),
        multivariate = FALSE
      )$psrf[, 1]
    }, error = function(e) {
      stats::setNames(rep(NA_real_, length(scalar_names)), scalar_names)
    })
  }

  list(
    rhats = rhats,
    ess = ess_vals,
    max_rhat = if (all(is.na(rhats))) NA_real_ else max(rhats, na.rm = TRUE),
    min_ess = if (all(is.na(ess_vals))) NA_real_ else min(ess_vals, na.rm = TRUE)
  )
}
