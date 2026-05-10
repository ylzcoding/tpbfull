library(MASS)

load_local_tpbfull <- function(path = ".") {
  env <- new.env(parent = globalenv())
  r_files <- list.files(file.path(path, "R"), pattern = "[.]R$", full.names = TRUE)
  for (f in r_files) {
    sys.source(f, envir = env)
  }
  env
}

sparse_data_gen <- function(n, p, num_nonzeros, seed,
                            signal_type = "uniform",
                            signal_params = list(df = 3, scale = 1),
                            target_snr = 4,
                            correlation = 0.9) {
  set.seed(seed)

  exponent <- abs(outer(1:p, 1:p, "-"))
  Sigma <- correlation^exponent
  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  X <- apply(X, 2, function(x) {(x - mean(x)) / sd(x)})

  nonzero_betas <- switch(
    signal_type,
    "uniform" = runif(num_nonzeros, signal_params$low, signal_params$up),
    "student_t" = rt(num_nonzeros, df = signal_params$df) * signal_params$scale,
    "bounded_uniform" = {
      mags <- runif(num_nonzeros, signal_params$low, signal_params$up)
      signs <- sample(c(-1, 1), num_nonzeros, replace = TRUE)
      mags * signs
    }
  )

  beta <- numeric(p)
  beta[sample(seq_len(p), size = num_nonzeros)] <- nonzero_betas

  signal_variance <- var(X %*% beta)
  sigmasq <- if (signal_variance < 1e-8) 1 else signal_variance / target_snr
  y <- c(X %*% beta + rnorm(n, 0, sqrt(sigmasq)))
  list(X = X, y = y, beta = beta, sigmasq = sigmasq)
}

summarize_scalars <- function(fit) {
  scalars <- fit$samples$scalars
  phi <- scalars[, "phi"]
  a <- scalars[, "a"]
  b <- scalars[, "b"]
  data.frame(
    phi_mean = mean(phi),
    phi_median = stats::median(phi),
    phi_q90 = unname(stats::quantile(phi, 0.90)),
    phi_q99 = unname(stats::quantile(phi, 0.99)),
    phi_max = max(phi),
    phi_gt_10 = mean(phi > 10),
    phi_gt_100 = mean(phi > 100),
    a_median = stats::median(a),
    b_median = stats::median(b),
    accept_a = fit$acceptance_rates$mean["a"],
    accept_b = fit$acceptance_rates$mean["b"],
    accept_phi = fit$acceptance_rates$mean["phi"],
    max_rhat = fit$diagnostics$max_rhat,
    min_ess = fit$diagnostics$min_ess
  )
}

run_phi_diagnostics <- function(seed = 1,
                                n = 50,
                                p = 200,
                                num_nonzeros = 4,
                                signal_type = "uniform",
                                signal_params = list(low = -5, up = 5),
                                target_snr = 4,
                                correlation = 0.9,
                                scale_phi_grid = c(0.05, 0.1, 0.5, 1),
                                configs = c("gamma_hs", "gamma_weak", "hcauchy_weak"),
                                num_chains = 2,
                                num_iter = 2000,
                                num_warmup = 1000,
                                thinning = 1,
                                proposal_type = "adaptive",
                                max_log_proposal_sd = 0.25,
                                woodbury = TRUE) {
  tpb <- load_local_tpbfull(".")
  data_gen <- sparse_data_gen(
    n = n, p = p, num_nonzeros = num_nonzeros, seed = seed,
    signal_type = signal_type, signal_params = signal_params,
    target_snr = target_snr, correlation = correlation
  )

  config_params <- list(
    gamma_hs = list(
      auto_find = FALSE,
      hyper_params = list(
        prior_type_a = "gamma", prior_type_b = "gamma",
        s_a = 1.5, r_a = 3,
        s_b = 1.5, r_b = 3
      )
    ),
    gamma_weak = list(
      auto_find = FALSE,
      hyper_params = list(
        prior_type_a = "gamma", prior_type_b = "gamma",
        s_a = 1.5, r_a = 0.25,
        s_b = 1.5, r_b = 0.25
      )
    ),
    hcauchy_weak = list(
      auto_find = FALSE,
      hyper_params = list(
        prior_type_a = "hcauchy", prior_type_b = "hcauchy",
        scale_a = 1, scale_b = 1
      )
    )
  )

  rows <- list()
  idx <- 1
  for (config_name in configs) {
    for (scale_phi in scale_phi_grid) {
      cfg <- config_params[[config_name]]
      cfg$hyper_params$scale_phi <- scale_phi
      set.seed(seed)
      fit <- tpb$tpb_full_pipeline(
        X = data_gen$X,
        y = data_gen$y,
        num_chains = num_chains,
        auto_find = cfg$auto_find,
        hyper_params = cfg$hyper_params,
        num_iter = num_iter,
        num_warmup = num_warmup,
        thinning = thinning,
        iter_pre_opt = 20,
        iter_selection = 100,
        woodbury = woodbury,
        proposal_type = proposal_type,
        max_log_proposal_sd = max_log_proposal_sd,
        model_prior_method = "complexity"
      )
      row <- summarize_scalars(fit)
      row$config <- config_name
      row$scale_phi <- scale_phi
      row$seed <- seed
      rows[[idx]] <- row
      idx <- idx + 1
      print(row)
    }
  }

  result <- do.call(rbind, rows)
  result <- result[, c("seed", "config", "scale_phi", setdiff(names(result), c("seed", "config", "scale_phi")))]
  print(result)
  invisible(result)
}

if (sys.nframe() == 0) {
  run_phi_diagnostics()
}
