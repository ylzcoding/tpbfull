library(MASS)

source("R/Gibbs_beta.R")
source("R/Gibbs_sigma.R")
source("R/Gibbs_zeta_psi.R")
source("R/log_marginal_beta.R")
source("R/mh_single.R")
source("R/fullGibbs.R")

sparse_data_gen <- function(n, p, num_nonzeros, seed,
                            signal_params = list(low = 0.5, up = 5),
                            target_snr = 4,
                            correlation = 0.8) {
  set.seed(seed)
  exponent <- abs(outer(seq_len(p), seq_len(p), "-"))
  Sigma <- correlation^exponent
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  X <- scale(X)

  mags <- runif(num_nonzeros, signal_params$low, signal_params$up)
  signs <- sample(c(-1, 1), num_nonzeros, replace = TRUE)
  beta <- numeric(p)
  beta[sample(seq_len(p), size = num_nonzeros)] <- mags * signs

  signal_variance <- var(as.vector(X %*% beta))
  sigmaSq <- if (signal_variance < 1e-8) 1 else signal_variance / target_snr
  y <- as.vector(X %*% beta + rnorm(n, sd = sqrt(sigmaSq)))
  list(X = X, y = y, beta = beta, sigmaSq = sigmaSq)
}

gamma_center <- function(a, b, shape = 1.5) {
  list(
    prior_type_a = "gamma",
    prior_type_b = "gamma",
    s_a = shape,
    r_a = shape / a,
    s_b = shape,
    r_b = shape / b,
    scale_a = 1,
    scale_b = 1,
    scale_phi = 1
  )
}

calc_metrics <- function(fit, beta_true) {
  scalars <- fit$samples$scalars
  beta_hat <- colMeans(fit$samples$beta)
  signal <- beta_true != 0
  phi <- scalars[, "phi"]
  a <- scalars[, "a"]
  b <- scalars[, "b"]

  data.frame(
    sse_total = sum((beta_true - beta_hat)^2),
    sse_signal = sum((beta_true[signal] - beta_hat[signal])^2),
    phi_mean = mean(phi),
    phi_median = median(phi),
    phi_q90 = as.numeric(quantile(phi, 0.90)),
    phi_q99 = as.numeric(quantile(phi, 0.99)),
    phi_gt_10 = mean(phi > 10),
    phi_gt_100 = mean(phi > 100),
    a_mean = mean(a),
    a_median = median(a),
    a_q90 = as.numeric(quantile(a, 0.90)),
    b_mean = mean(b),
    b_median = median(b),
    b_q90 = as.numeric(quantile(b, 0.90)),
    accept_a = fit$acceptance_rates$a,
    accept_b = fit$acceptance_rates$b,
    accept_phi = fit$acceptance_rates$phi,
    stringsAsFactors = FALSE
  )
}

run_one <- function(case, rho, seed, prior_name, prior,
                    proposal_type, mh_step_b = 0.2,
                    num_output = 1000, num_burnin = 1000) {
  data <- sparse_data_gen(
    n = case$n, p = case$p, num_nonzeros = case$k,
    seed = seed,
    target_snr = 4,
    correlation = rho
  )

  set.seed(seed)
  invisible(capture.output({
    fit <- fullGibbs(
      X = data$X,
      y = data$y,
      num_output = num_output,
      num_burnin = num_burnin,
      thin = 1,
      woodbury = TRUE,
      proposal_type = proposal_type,
      mh_step_b = mh_step_b,
      hyper_params = prior,
      max_log_proposal_sd = 0.5,
      max_log_mh_step = 1.0
    )
  }))

  cbind(
    data.frame(
      case = case$name,
      n = case$n,
      p = case$p,
      k = case$k,
      rho = rho,
      seed = seed,
      prior = prior_name,
      proposal_type = proposal_type,
      mh_step_b = mh_step_b,
      true_sigmaSq = data$sigmaSq,
      stringsAsFactors = FALSE
    ),
    calc_metrics(fit, data$beta)
  )
}

summarize_results <- function(rows) {
  aggregate(
    cbind(sse_total, sse_signal, phi_mean, phi_median, phi_q90, phi_q99,
          phi_gt_10, phi_gt_100, a_mean, a_median, a_q90,
          b_mean, b_median, b_q90, accept_a, accept_b, accept_phi) ~
      case + rho + prior + proposal_type + mh_step_b,
    data = rows,
    FUN = mean
  )
}

run_bi_tuning <- function(seeds = 1:2, mh_steps = c(0.1, 0.2, 0.3)) {
  case <- list(name = "mixed_p200_k40", n = 50, p = 200, k = 40)
  prior <- gamma_center(a = 0.5, b = 5)
  rows <- list()
  idx <- 1
  for (mh_step_b in mh_steps) {
    for (seed in seeds) {
      cat(sprintf("\nbi_adaptive tuning mh_step_b=%.2f seed=%d\n", mh_step_b, seed))
      rows[[idx]] <- run_one(
        case = case, rho = 0.4, seed = seed,
        prior_name = "normal_gamma_center",
        prior = prior,
        proposal_type = "bi_adaptive",
        mh_step_b = mh_step_b,
        num_output = 800,
        num_burnin = 800
      )
      print(rows[[idx]][, c("sse_signal", "phi_q90", "phi_q99",
                            "phi_gt_10", "phi_gt_100",
                            "b_q90", "accept_a", "accept_b", "accept_phi")])
      idx <- idx + 1
    }
  }
  result <- do.call(rbind, rows)
  write.csv(result, "bi_adaptive_b_step_tuning_results.csv", row.names = FALSE)
  summary <- summarize_results(result)
  write.csv(summary, "bi_adaptive_b_step_tuning_summary.csv", row.names = FALSE)
  summary
}

run_bi_vs_adaptive_grid <- function(seeds = 1:2,
                                    rhos = c(0.8, 0.4),
                                    bi_mh_step_b = 0.2) {
  cases <- list(
    sparse = list(name = "sparse_p200_k10", n = 50, p = 200, k = 10),
    mixed = list(name = "mixed_p200_k40", n = 50, p = 200, k = 40),
    dense = list(name = "dense_p100_k50", n = 50, p = 100, k = 50)
  )
  priors <- list(
    hs_center = gamma_center(a = 0.5, b = 0.5),
    normal_gamma_center = gamma_center(a = 0.5, b = 5),
    ridge_center = gamma_center(a = 5, b = 5)
  )

  rows <- list()
  idx <- 1
  for (case_name in names(cases)) {
    for (rho in rhos) {
      for (prior_name in names(priors)) {
        for (proposal_type in c("all_adaptive", "bi_adaptive")) {
          for (seed in seeds) {
            mh_step_b <- if (proposal_type == "bi_adaptive") bi_mh_step_b else 0.1
            cat(sprintf(
              "\ncase=%s rho=%.1f prior=%s proposal=%s seed=%d\n",
              cases[[case_name]]$name, rho, prior_name, proposal_type, seed
            ))
            rows[[idx]] <- run_one(
              case = cases[[case_name]],
              rho = rho,
              seed = seed,
              prior_name = prior_name,
              prior = priors[[prior_name]],
              proposal_type = proposal_type,
              mh_step_b = mh_step_b
            )
            print(rows[[idx]][, c("sse_signal", "phi_q90", "phi_q99",
                                  "phi_gt_10", "phi_gt_100",
                                  "a_q90", "b_q90",
                                  "accept_a", "accept_b", "accept_phi")])
            idx <- idx + 1
          }
        }
      }
    }
  }

  result <- do.call(rbind, rows)
  write.csv(result, "bi_vs_adaptive_results.csv", row.names = FALSE)
  summary <- summarize_results(result)
  write.csv(summary, "bi_vs_adaptive_summary.csv", row.names = FALSE)
  summary
}

if (sys.nframe() == 0) {
  print(run_bi_tuning())
  print(run_bi_vs_adaptive_grid())
}
