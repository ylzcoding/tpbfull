source("diagnostics/full_bayes_phi_diagnostics.R")

calculate_beta_metrics <- function(fit, beta_true) {
  beta_hat <- fit$samples$beta_summary$mean
  beta_sd <- fit$samples$beta_summary$sd
  lower <- fit$samples$beta_summary$lower
  upper <- fit$samples$beta_summary$upper
  is_signal <- beta_true != 0
  score <- abs(beta_hat) / (beta_sd + 1e-12)

  auc <- NA_real_
  if (any(is_signal) && any(!is_signal)) {
    wins <- outer(score[is_signal], score[!is_signal], ">")
    ties <- outer(score[is_signal], score[!is_signal], "==")
    auc <- (sum(wins) + 0.5 * sum(ties)) / (sum(is_signal) * sum(!is_signal))
  }

  c(
    sse = sum((beta_true - beta_hat)^2),
    auc = auc,
    csr = mean(((lower > 0) | (upper < 0)) == is_signal)
  )
}

summarize_fit <- function(fit, beta_true) {
  phi <- fit$scalar_traces$phi_samples
  a <- fit$scalar_traces$a_samples
  b <- fit$scalar_traces$b_samples
  c(
    calculate_beta_metrics(fit, beta_true),
    phi_median = stats::median(phi),
    phi_q90 = unname(stats::quantile(phi, 0.90)),
    phi_q99 = unname(stats::quantile(phi, 0.99)),
    phi_max = max(phi),
    phi_gt10 = mean(phi > 10),
    a_median = stats::median(a),
    b_median = stats::median(b),
    accept_phi = fit$acceptance_rates$mean["phi"]
  )
}

run_scale_sensitivity <- function(seed = 1,
                                  scale_phi_grid = c(0.05, 0.1, 0.5, 1),
                                  num_chains = 2,
                                  num_iter = 1200,
                                  num_warmup = 600,
                                  max_log_proposal_sd = 0.25) {
  tpb <- load_local_tpbfull(".")
  cases <- list(
    sparse = list(
      n = 50, p = 200, num_nonzeros = 4,
      signal_type = "uniform",
      signal_params = list(low = -5, up = 5),
      woodbury = TRUE,
      hyper_params = function(scale_phi) {
        list(
          prior_type_a = "gamma", prior_type_b = "gamma",
          s_a = 1.5, r_a = 3,
          s_b = 1.5, r_b = 3,
          scale_phi = scale_phi
        )
      }
    ),
    dense = list(
      n = 50, p = 50, num_nonzeros = 40,
      signal_type = "uniform",
      signal_params = list(low = 0, up = 1),
      woodbury = FALSE,
      hyper_params = function(scale_phi) {
        list(
          prior_type_a = "gamma", prior_type_b = "gamma",
          s_a = 1.5, r_a = 0.15,
          s_b = 1.5, r_b = 0.15,
          scale_phi = scale_phi
        )
      }
    )
  )

  rows <- list()
  idx <- 1
  for (case_name in names(cases)) {
    cfg <- cases[[case_name]]
    data_gen <- sparse_data_gen(
      n = cfg$n, p = cfg$p, num_nonzeros = cfg$num_nonzeros,
      seed = seed,
      signal_type = cfg$signal_type,
      signal_params = cfg$signal_params
    )

    for (scale_phi in scale_phi_grid) {
      set.seed(seed)
      fit <- tpb$tpb_full_pipeline(
        X = data_gen$X,
        y = data_gen$y,
        num_chains = num_chains,
        auto_find = FALSE,
        hyper_params = cfg$hyper_params(scale_phi),
        num_iter = num_iter,
        num_warmup = num_warmup,
        thinning = 1,
        woodbury = cfg$woodbury,
        proposal_type = "adaptive",
        max_log_proposal_sd = max_log_proposal_sd
      )
      row <- data.frame(
        case = case_name,
        seed = seed,
        scale_phi = scale_phi,
        t(summarize_fit(fit, data_gen$beta)),
        row.names = NULL
      )
      print(row)
      rows[[idx]] <- row
      idx <- idx + 1
    }
  }

  result <- do.call(rbind, rows)
  print(result)
  invisible(result)
}

if (sys.nframe() == 0) {
  run_scale_sensitivity()
}
