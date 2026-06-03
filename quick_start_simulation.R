library(tpbfull)
library(MASS)

simulate_sparse_regression <- function(n = 50, p = 200, s = 10,
                                       rho = 0.5, snr = 4, seed = 1) {
  set.seed(seed)

  Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  X <- scale(MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma))

  beta <- numeric(p)
  signal_id <- sample(seq_len(p), s)
  beta[signal_id] <- runif(s, -5, 5)

  signal <- as.vector(X %*% beta)
  sigmaSq <- stats::var(signal) / snr
  y <- signal + rnorm(n, sd = sqrt(sigmaSq))

  list(X = X, y = y, beta = beta, sigmaSq = sigmaSq)
}

beta_metrics <- function(beta_summary, beta_true) {
  beta_hat <- beta_summary$mean
  selected <- beta_summary$lower > 0 | beta_summary$upper < 0
  true_signal <- beta_true != 0

  c(
    mse = mean((beta_hat - beta_true)^2),
    sse_signal = sum((beta_hat[true_signal] - beta_true[true_signal])^2),
    selection_accuracy = mean(selected == true_signal)
  )
}

dat <- simulate_sparse_regression(
  n = 50,
  p = 200,
  s = 10,
  rho = 0.5,
  snr = 9,
  seed = 1
)

fit <- tpb_full_pipeline(
  X = dat$X,
  y = dat$y,
  num_chains = 4,
  num_iter = 4000,
  num_warmup = 2000,
  thinning = 1,
  auto_find = TRUE,
  init_method = "ridge",
  iter_pre_opt = 20,
  pre_opt_burnin = 50,
  pre_opt_samples = 50,
  selection_samples = 100,
  woodbury = TRUE,
  verbose = FALSE
)

cat("Winner:", fit$auto_find$winner_name, "\n")
print(fit$auto_find$winning_modes)
print(round(fit$auto_find$model_probabilities, 3))
print(round(beta_metrics(fit$samples$beta_summary, dat$beta), 4))
