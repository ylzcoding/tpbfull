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
    "uniform" = {
      runif(num_nonzeros, signal_params$low, signal_params$up)
    },
    "student_t" = {
      rt(num_nonzeros, df = signal_params$df) * signal_params$scale
    },
    "bounded_uniform" = {
      mags <- runif(num_nonzeros, signal_params$low, signal_params$up)
      signs <- sample(c(-1, 1), num_nonzeros, replace = TRUE)
      mags * signs
    }
  )

  beta <- numeric(p)
  signal_indices <- sample(1:p, size = num_nonzeros)
  beta[signal_indices] <- nonzero_betas

  signal_variance <- var(X %*% beta)
  sigmasq <- if (signal_variance < 1e-8) {
    1.0
  } else {
    signal_variance / target_snr
  }

  y <- c(X %*% beta + rnorm(n, 0, sqrt(sigmasq)))
  list(X = X, y = y, beta = beta, sigmasq = sigmasq)
}

run_repro <- function(seeds = 1:15,
                      n = 50,
                      p = 100,
                      num_nonzeros = 3,
                      signal_type = "uniform",
                      signal_params = list(low = 0, up = 1),
                      correlation = 0.9,
                      target_snr = 4,
                      model_prior_method = "complexity",
                      model_size_prior_alpha = 1,
                      model_size_prior_beta = 1,
                      candidates = list(
                        hs = list(a = 0.5, b = 0.5),
                        normal_gamma = list(a = 1, b = 10),
                        ridge = list(a = 10, b = 10)
                      ),
                      iter_pre_opt = 50,
                      iter_selection = 1000,
                      pre_opt_burnin = 200,
                      pre_opt_samples = 200,
                      woodbury = TRUE) {
  tpb <- load_local_tpbfull(".")
  rows <- vector("list", length(seeds))

  for (i in seq_along(seeds)) {
    seed <- seeds[i]
    data_gen <- sparse_data_gen(
      n = n, p = p, num_nonzeros = num_nonzeros,
      seed = seed,
      signal_type = signal_type,
      signal_params = signal_params,
      target_snr = target_snr,
      correlation = correlation
    )

    set.seed(seed)
    out <- capture.output({
      fit <- tpb$run_model_competition(
        X = data_gen$X,
        y = data_gen$y,
        iter_pre_opt = iter_pre_opt,
        iter_selection = iter_selection,
        pre_opt_burnin = pre_opt_burnin,
        pre_opt_samples = pre_opt_samples,
        woodbury = woodbury,
        candidates = candidates,
        use_model_prior = TRUE,
        model_prior_method = model_prior_method,
        model_size_prior_alpha = model_size_prior_alpha,
        model_size_prior_beta = model_size_prior_beta
      )
    })

    probs <- round(fit$raw$model_probabilities, 4)
    log_prior <- fit$model_prior$log_prior
    rows[[i]] <- c(
      seed = seed,
      active_prop = fit$model_prior$active_prop,
      active_count = fit$model_prior$active_count,
      winner = fit$winner_name,
      stats::setNames(as.character(probs), paste0(names(probs), "_prob")),
      stats::setNames(as.character(log_prior), paste0(names(log_prior), "_log_prior"))
    )
  }

  result <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  result$seed <- as.integer(result$seed)
  result$active_prop <- as.numeric(result$active_prop)
  numeric_cols <- setdiff(names(result), c("seed", "winner"))
  result[numeric_cols] <- lapply(result[numeric_cols], as.numeric)
  print(result)
  print(table(result$winner))
  invisible(result)
}

if (sys.nframe() == 0) {
  run_repro()
}
