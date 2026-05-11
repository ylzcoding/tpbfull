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
                            signal_params = list(low = 0, up = 1),
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
    "bounded_uniform" = {
      mags <- runif(num_nonzeros, signal_params$low, signal_params$up)
      signs <- sample(c(-1, 1), num_nonzeros, replace = TRUE)
      mags * signs
    },
    "student_t" = rt(num_nonzeros, df = signal_params$df) * signal_params$scale
  )

  beta <- numeric(p)
  beta[sample.int(p, size = num_nonzeros)] <- nonzero_betas
  signal_variance <- var(X %*% beta)
  sigmasq <- if (signal_variance < 1e-8) 1 else signal_variance / target_snr
  y <- c(X %*% beta + rnorm(n, 0, sqrt(sigmasq)))
  list(X = X, y = y, beta = beta, sigmasq = sigmasq)
}

candidate_sets <- list(
  ridge5 = list(
    hs = list(a = 0.5, b = 0.5),
    normal_gamma = list(a = 1, b = 10),
    ridge = list(a = 5, b = 5)
  ),
  ridge = list(
    hs = list(a = 0.5, b = 0.5),
    normal_gamma = list(a = 1, b = 10),
    ridge = list(a = 10, b = 10)
  ),
  ridge20 = list(
    hs = list(a = 0.5, b = 0.5),
    normal_gamma = list(a = 1, b = 10),
    ridge = list(a = 20, b = 20)
  ),
  ridge_student_t = list(
    hs = list(a = 0.5, b = 0.5),
    normal_gamma = list(a = 1, b = 10),
    ridge = list(a = 10, b = 10),
    student_t = list(a = 10, b = 1)
  ),
  student_t = list(
    hs = list(a = 0.5, b = 0.5),
    normal_gamma = list(a = 1, b = 10),
    student_t = list(a = 10, b = 1)
  ),
  sparse_ng_student_t = list(
    hs = list(a = 0.5, b = 0.5),
    normal_gamma_sparse = list(a = 0.5, b = 10),
    student_t = list(a = 10, b = 1)
  ),
  four_way = list(
    hs = list(a = 0.5, b = 0.5),
    normal_gamma_sparse = list(a = 0.5, b = 10),
    normal_gamma = list(a = 1, b = 10),
    student_t = list(a = 10, b = 1)
  )
)

cases <- list(
  dense_50_50_40_weak = list(
    n = 50, p = 50, k = 40, signal_type = "uniform",
    signal_params = list(low = 0, up = 1), target_snr = 4, correlation = 0.9,
    target = "dense"
  ),
  sparse_50_200_10_weak = list(
    n = 50, p = 200, k = 10, signal_type = "uniform",
    signal_params = list(low = 0, up = 1), target_snr = 4, correlation = 0.9,
    target = "sparse"
  ),
  dense_50_200_60_weak = list(
    n = 50, p = 200, k = 60, signal_type = "uniform",
    signal_params = list(low = 0, up = 1), target_snr = 4, correlation = 0.9,
    target = "dense"
  ),
  sparse_50_200_10_strong = list(
    n = 50, p = 200, k = 10, signal_type = "bounded_uniform",
    signal_params = list(low = 0.5, up = 5), target_snr = 4, correlation = 0.9,
    target = "sparse"
  ),
  dense_50_200_60_strong = list(
    n = 50, p = 200, k = 60, signal_type = "bounded_uniform",
    signal_params = list(low = 0.5, up = 5), target_snr = 4, correlation = 0.9,
    target = "dense"
  )
)

winner_family <- function(winner, candidates) {
  ab <- candidates[[winner]]
  if (is.null(ab)) {
    return(NA_character_)
  }
  if (ab$a <= 0.5) {
    "sparse"
  } else {
    "dense"
  }
}

run_one_setting <- function(tpb, case, candidates, seeds,
                            threshold, alpha, beta,
                            selection_score,
                            iter_pre_opt, iter_selection,
                            pre_opt_burnin, pre_opt_samples) {
  rows <- vector("list", length(seeds))
  for (i in seq_along(seeds)) {
    seed <- seeds[i]
    data <- sparse_data_gen(
      n = case$n, p = case$p, num_nonzeros = case$k, seed = seed,
      signal_type = case$signal_type,
      signal_params = case$signal_params,
      target_snr = case$target_snr,
      correlation = case$correlation
    )
    set.seed(seed)
    invisible(capture.output({
      fit <- tpb$run_model_competition(
        X = data$X, y = data$y,
        iter_pre_opt = iter_pre_opt,
        iter_selection = iter_selection,
        pre_opt_burnin = pre_opt_burnin,
        pre_opt_samples = pre_opt_samples,
        woodbury = TRUE,
        candidates = candidates,
        use_model_prior = TRUE,
        sparsity_threshold = threshold,
        selection_score = selection_score,
        model_prior_method = "binomial",
        model_size_prior_alpha = alpha,
        model_size_prior_beta = beta
      )
    }))
    rows[[i]] <- data.frame(
      seed = seed,
      active_count = fit$model_prior$active_count,
      active_prop = fit$model_prior$active_prop,
      winner = fit$winner_name,
      winner_family = winner_family(fit$winner_name, candidates),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_setting <- function(rows, target) {
  winners <- paste(names(table(rows$winner)), as.integer(table(rows$winner)),
                   sep = ":", collapse = ";")
  data.frame(
    target = target,
    n = nrow(rows),
    correct_rate = mean(rows$winner_family == target),
    active_mean = mean(rows$active_count),
    active_min = min(rows$active_count),
    active_max = max(rows$active_count),
    winners = winners,
    stringsAsFactors = FALSE
  )
}

run_grid <- function(
    seeds = 1:5,
    case_names = names(cases),
    candidate_names = names(candidate_sets),
    thresholds = c(0.05, 0.075, 0.10, 0.15),
    beta_priors = list(c(0.5, 0.5), c(1, 1), c(0.5, 2), c(2, 0.5)),
    selection_scores = c("data_likelihood"),
    iter_pre_opt = 15,
    iter_selection = 200,
    pre_opt_burnin = 40,
    pre_opt_samples = 40,
    out_csv = NULL) {
  tpb <- load_local_tpbfull(".")
  summaries <- list()
  idx <- 1
  for (case_name in case_names) {
    case <- cases[[case_name]]
    for (candidate_name in candidate_names) {
      candidates <- candidate_sets[[candidate_name]]
      for (threshold in thresholds) {
        for (prior in beta_priors) {
          for (selection_score in selection_scores) {
          alpha <- prior[[1]]
          beta <- prior[[2]]
          cat(sprintf(
            "\ncase=%s candidates=%s threshold=%.3f alpha=%.2f beta=%.2f score=%s\n",
            case_name, candidate_name, threshold,
            alpha, beta, selection_score
          ))
          rows <- run_one_setting(
            tpb = tpb, case = case, candidates = candidates, seeds = seeds,
            threshold = threshold, alpha = alpha, beta = beta,
            selection_score = selection_score,
            iter_pre_opt = iter_pre_opt, iter_selection = iter_selection,
            pre_opt_burnin = pre_opt_burnin, pre_opt_samples = pre_opt_samples
          )
          smry <- summarize_setting(rows, case$target)
          smry$case <- case_name
          smry$candidate_set <- candidate_name
          smry$threshold <- threshold
          smry$alpha <- alpha
          smry$beta <- beta
          smry$selection_score <- selection_score
          summaries[[idx]] <- smry
          idx <- idx + 1
          print(smry)
          }
        }
      }
    }
  }
  result <- do.call(rbind, summaries)
  result <- result[order(-result$correct_rate, result$case,
                         result$candidate_set, result$threshold), ]
  if (!is.null(out_csv)) {
    utils::write.csv(result, out_csv, row.names = FALSE)
  }
  result
}

if (sys.nframe() == 0) {
  result <- run_grid(out_csv = "diagnostics/binomial_prior_grid_results.csv")
  print(head(result, 30))
}
