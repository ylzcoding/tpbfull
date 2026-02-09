library(tpbfull)
library(MASS)
library(rstan)
# devtools::install_github("ylzcoding/tpbfull", ref = "main", upgrade = "never")

# optimize rstan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

calculate_metrics <- function(beta_samples, beta_true) {
  # beta_samples dimension: num_samples * p
  beta_hat <- colMeans(beta_samples)
  beta_sd <- apply(beta_samples, 2, sd)

  # Calculate total SSE
  sse_total <- sum((beta_true - beta_hat)^2)

  # Calculate SSE for non-zero coefficients (signals)
  signal_idx <- which(beta_true != 0)
  noise_idx <- which(beta_true == 0)

  sse_signal <- if (length(signal_idx) > 0) {
    sum((beta_true[signal_idx] - beta_hat[signal_idx])^2)
  } else {
    0 # No signals, so no signal error
  }

  # Calculate average 95% credible interval width and FPR for zero coefficients (noise)
  ci_bounds <- apply(beta_samples, 2, quantile, probs = c(0.025, 0.975))
  lower <- ci_bounds[1, ]
  upper <- ci_bounds[2, ]

  is_selected <- (lower > 0) | (upper < 0) # "Selected" if 0 is NOT in the interval (estimated as signals)
  is_true_signal <- (beta_true != 0)

  # csr, correct selection rate = Sum of correct decisions / Total number of coefficients
  correct_selection <- (is_selected == is_true_signal)
  csr <- mean(correct_selection)

  # auc
  predictor_score <- abs(beta_hat) / (beta_sd + 1e-12)
  true_labels <- as.numeric(is_true_signal)

  scores_signal <- predictor_score[true_labels == 1]
  scores_noise  <- predictor_score[true_labels == 0]
  n_signal <- length(scores_signal)
  n_noise  <- length(scores_noise)
  wins <- outer(scores_signal, scores_noise, ">")
  ties <- outer(scores_signal, scores_noise, "==")
  auc_val <- (sum(wins) + 0.5 * sum(ties)) / (n_signal * n_noise)

  return(list(
    sse_total = sse_total,
    sse_signal = sse_signal,
    auc = auc_val,
    csr = csr
  ))
}



rep_func <- function(X, y, beta_true, stan_file_path, num_chains, num_iter, num_warmup, thinning, woodbury, seed){
  # number of samples we'd like to collect
  # X should a design matrix
  # y should be a numeric vector

  num_output <- num_iter - num_warmup
  num_burnin <- num_warmup
  set.seed(seed)

  # =========================================================================
  # RStan (Half-Cauchy Prior)
  # =========================================================================
  cat("Running Stan (Gamma Priors)...\n")

  stan_data <- list(
    N = nrow(X),
    p = ncol(X),
    X = X,
    y = as.vector(y)
  )

  fit_stan <- stan(
    file = stan_file_path,
    data = stan_data,
    chains = num_chains,
    iter = num_iter,
    warmup = num_warmup,
    seed = seed,
    cores = min(num_chains, parallel::detectCores()),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    refresh = 1000
  )

  # Extract samples
  extracted <- rstan::extract(fit_stan)
  beta_samples_stan <- extracted$beta
  stan_metrics <- calculate_metrics(beta_samples_stan, beta_true)

  hyper_pars_stan <- data.frame(
    a = extracted$a,
    b = extracted$b,
    phi = extracted$phi,
    sigmaSq = extracted$sigmaSq
  )

  fit_summary_stan <- rstan::summary(fit_stan, pars = c("a", "b", "phi", "sigmaSq"), probs = c())$summary
  stan_results <- list(
    a_samples = as.vector(hyper_pars_stan[,"a"]),
    b_samples = as.vector(hyper_pars_stan[,"b"]),
    phi_samples = as.vector(hyper_pars_stan[,"phi"]),
    metrics = stan_metrics,
    hyper_means = colMeans(hyper_pars_stan),
    convergence = list(
      max_rhat = max(fit_summary_stan[, "Rhat"], na.rm = TRUE),
      min_ess  = min(fit_summary_stan[, "n_eff"], na.rm = TRUE),
      divergences = rstan::get_num_divergent(fit_stan)
    )
  )


  # =========================================================================
  # Gibbs Sampler (Gamma Prior on a, b)
  # =========================================================================
  cat("Running Gibbs Sampler (Gamma Priors)...\n")
  # Run multiple chains to calculate R-hat
  gibbs_chains_scalars <- list() # Store scalar samples for diagnostics
  gibbs_chains_beta    <- list() # Store beta samples
  for (i in 1:num_chains) {
    cat(sprintf("  Chain %d / %d\n", i, num_chains))
    chain_res <- fullGibbs(
      X = X, y = y,
      num_output = num_output,
      num_burnin = num_burnin,
      thin = thinning,
      woodbury = woodbury,
      hyper_params = list(s_a=2, r_a=2, s_b=2, r_b=1)
    )
    gibbs_chains_scalars[[i]] <- chain_res$samples$scalars
    gibbs_chains_beta[[i]]    <- chain_res$samples$beta
  }

  # R-hat
  mcmc_list_scalars <- coda::mcmc.list(lapply(gibbs_chains_scalars, coda::mcmc))
  gelman_diag <- tryCatch({
    coda::gelman.diag(mcmc_list_scalars, multivariate = FALSE)$psrf[, 1]
  }, error = function(e) return(rep(NA, ncol(gibbs_chains_scalars[[1]]))))

  # Calculate ESS (Effective Sample Size)
  ess_vals <- tryCatch({
    coda::effectiveSize(mcmc_list_scalars)
  }, error = function(e) return(rep(NA, ncol(gibbs_chains_scalars[[1]]))))

  # Stack all chains into one matrix: (num_chains * num_output) x p
  gibbs_beta_combined <- do.call(rbind, gibbs_chains_beta)
  gibbs_scalar_combined <- do.call(rbind, gibbs_chains_scalars)
  gibbs_metrics <- calculate_metrics(gibbs_beta_combined, beta_true)
  gibbs_hyper_means <- colMeans(gibbs_scalar_combined)

  gibbs_results <- list(
    a_samples = as.vector(gibbs_scalar_combined[,"a"]),
    b_samples = as.vector(gibbs_scalar_combined[,"b"]),
    phi_samples = as.vector(gibbs_scalar_combined[,"phi"]),
    metrics = gibbs_metrics,
    hyper_means = gibbs_hyper_means,
    convergence = list(
      max_rhat = max(gelman_diag, na.rm = TRUE),
      min_ess  = min(ess_vals, na.rm = TRUE),
      divergences = NA, # Gibbs does not produce HMC divergences
      rhats = gelman_diag # Detailed R-hats
    )
  )

  return(list(
    gamma_gibbs = gibbs_results,
    gamma_stan = stan_results
  ))
}

sparse_data_gen <- function(n, p, num_nonzeros, seed,
                            signal_type = "uniform", # student_t of uniform
                            signal_params = list(df = 3, scale = 1),
                            target_snr = 4,
                            correlation = 0.9) {

  set.seed(seed)

  # Sigma[i, j] = correlation^|i-j|
  exponent <- abs(outer(1:p, 1:p, "-"))
  Sigma <- correlation^exponent

  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  X <- apply(X, 2, function(x) {(x - mean(x)) / sd(x)})

  nonzero_betas <- switch(
    signal_type,
    "uniform" = {
      runif(num_nonzeros, signal_params$low, signal_params$up)
    },
    "student_t" = {
      rt(num_nonzeros, df = signal_params$df) * signal_params$scale
    },
    runif(num_nonzeros, -5, 5)
  )

  beta <- numeric(p)
  signal_indices <- sample(1:p, size = num_nonzeros)
  beta[signal_indices] <- nonzero_betas

  # keep SNR fixed
  signal_variance <- var(X %*% beta)

  if (signal_variance < 1e-8) {
    sigmasq <- 1.0
  } else {
    sigmasq <- signal_variance / target_snr # min(max(1.0, signal_variance / target_snr), 100.0)
  }

  y <- c(X %*% beta + rnorm(n, 0, sqrt(sigmasq)))
  return(list(X = X, y = y, beta = beta, sigmasq = sigmasq))
}

######### SLURM settings ###########

dir <- "my_dir"
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#cat("Running seed ", task_id, "\n")
n <- 50
p <- 200
num_nonzeros <- 4 # 4, 16, 30, or 60
unif_low <- -5 # -1 & -5
unif_up <- 5 # 1 & 5
df <- 3
scale <- 1
corr <- 0.9 # low 0.5 & super 0.9
target_snr <- 4


data_gen <- sparse_data_gen(
  n = n, p = p,
  num_nonzeros = num_nonzeros,
  seed = task_id,
  signal_type = "uniform",
  signal_params = list(low = unif_low, up = unif_up),   #list(df = df, scale = scale),list(low = unif_low, up = unif_up)
  target_snr = target_snr, # or 9
  correlation = corr
)


X <- data_gen$X
y <- data_gen$y
beta_true <- data_gen$beta

sim_res <- rep_func(X = X, y = y, beta_true = beta_true,
                    stan_file_path = "tpb_with_gamma.stan",
                    num_chains = 5, num_iter = 2000, num_warmup = 1000, thinning = 1,
                    woodbury = TRUE, seed = task_id)
