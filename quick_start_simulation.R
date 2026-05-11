library(tpbfull)
library(MASS)

# This script gives a small simulation example for tpb_full_pipeline().

sparse_data_gen <- function(n, p, num_nonzeros, seed,
                            signal_type = "bounded_uniform",
                            signal_params = list(low = 0.5, up = 5),
                            target_snr = 4,
                            correlation = 0.8) {
  set.seed(seed)

  exponent <- abs(outer(seq_len(p), seq_len(p), "-"))
  Sigma <- correlation^exponent

  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  X <- scale(X)

  nonzero_betas <- switch(
    signal_type,
    uniform = {
      runif(num_nonzeros, signal_params$low, signal_params$up)
    },
    bounded_uniform = {
      mags <- runif(num_nonzeros, signal_params$low, signal_params$up)
      signs <- sample(c(-1, 1), num_nonzeros, replace = TRUE)
      mags * signs
    },
    student_t = {
      rt(num_nonzeros, df = signal_params$df) * signal_params$scale
    },
    stop("Unknown signal_type.")
  )

  beta <- numeric(p)
  signal_idx <- sample(seq_len(p), size = num_nonzeros)
  beta[signal_idx] <- nonzero_betas

  signal_variance <- var(as.vector(X %*% beta))
  sigmaSq <- if (signal_variance < 1e-8) {
    1
  } else {
    signal_variance / target_snr
  }

  y <- as.vector(X %*% beta + rnorm(n, sd = sqrt(sigmaSq)))
  list(X = X, y = y, beta = beta, sigmaSq = sigmaSq)
}

calculate_metrics <- function(beta_summary, beta_true) {
  beta_hat <- beta_summary$mean
  beta_sd <- beta_summary$sd
  selected <- beta_summary$lower > 0 | beta_summary$upper < 0
  true_signal <- beta_true != 0

  list(
    sse_total = sum((beta_true - beta_hat)^2),
    sse_signal = sum((beta_true[true_signal] - beta_hat[true_signal])^2),
    csr = mean(selected == true_signal),
    signal_rank_score = abs(beta_hat) / (beta_sd + 1e-12)
  )
}

# Example cases:
#   sparse high-dimensional: n = 50, p = 200, num_nonzeros = 10
#   dense lower-dimensional: n = 50, p = 50,  num_nonzeros = 35
data <- sparse_data_gen(
  n = 50,
  p = 200,
  num_nonzeros = 10,
  seed = 1,
  signal_type = "bounded_uniform",
  signal_params = list(low = 0.5, up = 5),
  target_snr = 4,
  correlation = 0.8
)

# Model-selection scores:
#   "beta_prior":
#      original score, log p(beta | a, b, phi)
#   "posterior_kernel":
#      unnormalized beta posterior kernel,
#      log p(y | X, beta, sigmaSq) + log p(beta | a, b, phi)
#
# The empirical-Bayes competition first selects a candidate (a,b), then the
# fully Bayesian sampler uses Gamma priors centered at the winning values.
fit <- tpb_full_pipeline(
  X = data$X,
  y = data$y,
  num_chains = 4,
  num_iter = 4000,
  num_warmup = 2000,
  thinning = 1,
  auto_find = TRUE,
  selection_score = "posterior_kernel",
  iter_pre_opt = 20,
  pre_opt_burnin = 50,
  pre_opt_samples = 50,
  iter_selection = 500,
  woodbury = TRUE
)

fit$auto_find$winner_name
fit$auto_find$winning_modes
fit$auto_find$model_probabilities

metrics <- calculate_metrics(fit$samples$beta_summary, data$beta)
metrics[c("sse_total", "sse_signal", "csr")]
