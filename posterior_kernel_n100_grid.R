library(MASS)

source("R/Gibbs_beta.R")
source("R/Gibbs_zeta_psi.R")
source("R/model_competition.R")

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
    uniform = runif(num_nonzeros, signal_params$low, signal_params$up),
    bounded_uniform = {
      mags <- runif(num_nonzeros, signal_params$low, signal_params$up)
      signs <- sample(c(-1, 1), num_nonzeros, replace = TRUE)
      mags * signs
    },
    student_t = rt(num_nonzeros, df = signal_params$df) * signal_params$scale,
    stop("Unknown signal_type.")
  )

  beta <- numeric(p)
  beta[sample(seq_len(p), size = num_nonzeros)] <- nonzero_betas
  signal_variance <- var(as.vector(X %*% beta))
  sigmaSq <- if (signal_variance < 1e-8) 1 else signal_variance / target_snr
  y <- as.vector(X %*% beta + rnorm(n, sd = sqrt(sigmaSq)))

  list(X = X, y = y, beta = beta, sigmaSq = sigmaSq)
}

winner_family <- function(winner) {
  if (winner == "hs") {
    "sparse"
  } else {
    "dense"
  }
}

run_one_case <- function(n, p, k, rho, candidates, seeds,
                         iter_pre_opt = 8,
                         iter_selection = 100,
                         pre_opt_burnin = 30,
                         pre_opt_samples = 30) {
  rows <- vector("list", length(seeds))
  for (i in seq_along(seeds)) {
    seed <- seeds[i]
    data <- sparse_data_gen(
      n = n, p = p, num_nonzeros = k, seed = seed,
      signal_type = "bounded_uniform",
      signal_params = list(low = 0.5, up = 5),
      target_snr = 4,
      correlation = rho
    )

    set.seed(seed)
    fit <- suppressMessages(suppressWarnings(capture.output({
      out <- run_model_competition(
        X = data$X,
        y = data$y,
        iter_pre_opt = iter_pre_opt,
        iter_selection = iter_selection,
        pre_opt_burnin = pre_opt_burnin,
        pre_opt_samples = pre_opt_samples,
        woodbury = TRUE,
        candidates = candidates,
        selection_score = "posterior_kernel"
      )
    })))

    rows[[i]] <- data.frame(
      seed = seed,
      winner = out$winner_name,
      family = winner_family(out$winner_name),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_case <- function(rows, target) {
  tab <- table(rows$winner)
  data.frame(
    target = target,
    n_rep = nrow(rows),
    sparse_rate = mean(rows$family == "sparse"),
    dense_rate = mean(rows$family == "dense"),
    target_rate = if (target == "mixed") NA_real_ else mean(rows$family == target),
    winners = paste(names(tab), as.integer(tab), sep = ":", collapse = ";"),
    stringsAsFactors = FALSE
  )
}

run_grid <- function(seeds = 1:10,
                     correlations = c(0.8, 0.4),
                     ridge_values = c(5, 10)) {
  cases <- list(
    p200_k10 = list(n = 100, p = 200, k = 10, target = "sparse"),
    p200_k20 = list(n = 100, p = 200, k = 20, target = "sparse"),
    p200_k40 = list(n = 100, p = 200, k = 40, target = "mixed"),
    p100_k30 = list(n = 100, p = 100, k = 30, target = "mixed"),
    p100_k50 = list(n = 100, p = 100, k = 50, target = "dense"),
    p100_k70 = list(n = 100, p = 100, k = 70, target = "dense")
  )

  out <- list()
  idx <- 1
  for (ridge in ridge_values) {
    candidates <- list(
      hs = list(a = 0.5, b = 0.5),
      normal_gamma = list(a = 1, b = 10),
      ridge = list(a = ridge, b = ridge)
    )
    for (rho in correlations) {
      for (case_name in names(cases)) {
        case <- cases[[case_name]]
        cat(sprintf(
          "\ncase=%s rho=%.1f ridge=%s\n", case_name, rho, ridge
        ))
        rows <- run_one_case(
          n = case$n, p = case$p, k = case$k, rho = rho,
          candidates = candidates, seeds = seeds
        )
        smry <- summarize_case(rows, case$target)
        smry$case <- case_name
        smry$rho <- rho
        smry$ridge <- ridge
        out[[idx]] <- smry
        idx <- idx + 1
        print(smry)
      }
    }
  }

  result <- do.call(rbind, out)
  result <- result[order(result$ridge, result$rho, result$case), ]
  write.csv(result, "posterior_kernel_n100_grid_results.csv", row.names = FALSE)
  result
}

if (sys.nframe() == 0) {
  print(run_grid())
}
