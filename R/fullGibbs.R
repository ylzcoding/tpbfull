#' Fully Bayesian TPB sampler with marginal and collapsed updates
#'
#' @import mvtnorm
#' @import coda
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param num_output Integer, number of post-warmup iterations before thinning
#' @param num_burnin Integer, number of burn-in iterations to discard
#' @param thin Integer, thinning interval for saved posterior samples
#' @param woodbury Logical, use Woodbury identity in beta update
#' @param diagX Logical, assume diagonal X
#' @param adapt_block_size Integer, number of iterations per adaptation block
#' @param r_opt Numeric, target acceptance rate for adaptive MH (default 0.3)
#' @param max_log_proposal_sd Maximum marginal proposal standard deviation on
#'   the log scale for adaptive MH covariance regularization.
#' @param max_log_mh_step Maximum single proposal move on the log scale. Proposals
#'   beyond this distance are treated as self-loops.
#' @param sampler "marginal" integrates local scales out of the hyperparameter
#'   update; "conditional" retains the original nu, lambda and xi augmentation.
#' @param interweave_every Apply the additional block once every this many iterations.
#' @param adapt_slice Tune slice widths during burn-in, then freeze them.
#' @param interweave Add a beta-collapsed shape/scale update in transformed local scales.
#' @param hyper_params List of hyperparameters (prior_type_a, prior_type_b,
#'   s_a, r_a, s_b, r_b, scale_a, scale_b, prior_type_phi, scale_phi).
#'   prior_type_phi = "sqrt_hcauchy" (default) puts a half-Cauchy prior on
#'   sqrt(phi); "hcauchy" retains the earlier prior on phi itself. scale_phi
#'   is the half-Cauchy scale of the selected quantity. If scale_phi is NULL
#'   or invalid, scale_phi = 1 is used.
#' @details The marginal sampler jointly proposes log(a) and log(phi), then
#'   updates log(b) with a random-walk proposal of standard deviation 0.1.
#'   Its proposal covariance is adapted during burn-in and fixed afterward.
#'   The additional interweaving block updates a, b and phi after integrating
#'   beta out, then immediately redraws beta. It preserves the model and priors.
#'   Set interweave = FALSE to skip this block; local scales still use log storage.
#' @return Posterior samples, MH acceptance rates and final tuning values.
#'   beta_loglik is log p(beta | a, b, phi); total_logpost adds the hyperpriors.
#'   For the conditional sampler, marginal log densities are NULL and MH
#'   acceptance rates are NA because its hyperparameters use slice updates.
#' @export
fullGibbs <- function(X, y, num_output = 10000, num_burnin = 10000, thin = 1,
                      woodbury = TRUE, diagX = FALSE,
                      adapt_block_size = 100, r_opt = 0.3,
                      max_log_proposal_sd = 0.5,
                      max_log_mh_step = 1.0,
                      hyper_params = list(prior_type_a = "gamma", prior_type_b = "gamma",
                                          s_a = 2.5, r_a = 5, s_b = 2.5, r_b = 2.5,
                                          scale_a = 1, scale_b = 1,
                                          prior_type_phi = "sqrt_hcauchy", scale_phi = 1),
                      interweave = sampler[1] == "marginal",
                      sampler = c("marginal", "conditional"),
                      interweave_every = 1,
                      adapt_slice = FALSE) {

  sampler <- match.arg(sampler)
  if (sampler == "conditional" && interweave)
    stop("Use interweave = FALSE with the conditional sampler.")
  if (length(interweave_every) != 1 || !is.finite(interweave_every) ||
      interweave_every < 1 || interweave_every != floor(interweave_every))
    stop("interweave_every must be a positive integer.")
  slice_width <- rep(1, 3)
  slice_jumps <- numeric(3)
  slice_updates <- 0L
  n <- nrow(X)
  p <- ncol(X)
  total_iter <- num_output + num_burnin
  n_save <- floor(num_output / thin)
  hyper_params <- modifyList(
    list(prior_type_a = "gamma", prior_type_b = "gamma",
         s_a = 2.5, r_a = 5, s_b = 2.5, r_b = 2.5,
         scale_a = 1, scale_b = 1,
         prior_type_phi = "sqrt_hcauchy", scale_phi = 1),
    hyper_params
  )
  hyper_params$prior_type_a <- match.arg(hyper_params$prior_type_a, c("gamma", "hcauchy"))
  hyper_params$prior_type_b <- match.arg(hyper_params$prior_type_b, c("gamma", "hcauchy"))
  hyper_params$prior_type_phi <- match.arg(hyper_params$prior_type_phi,
                                          c("sqrt_hcauchy", "hcauchy"))
  if (is.null(hyper_params$scale_phi) ||
      !is.finite(hyper_params$scale_phi) ||
      hyper_params$scale_phi <= 0) {
    hyper_params$scale_phi <- 1
  }

  regularize_log_proposal_cov <- function(cov_mat) {
    cov_mat <- as.matrix(cov_mat)
    cov_mat <- (cov_mat + t(cov_mat)) / 2
    diag(cov_mat) <- pmax(diag(cov_mat), 1e-8)

    current_max_sd <- sqrt(max(diag(cov_mat)))
    if (is.finite(current_max_sd) && current_max_sd > max_log_proposal_sd) {
      cov_mat <- cov_mat * (max_log_proposal_sd / current_max_sd)^2
      diag(cov_mat) <- pmax(diag(cov_mat), 1e-8)
    }

    eig <- eigen(cov_mat, symmetric = TRUE)
    if (min(eig$values) <= 1e-8) {
      cov_mat <- eig$vectors %*% diag(pmax(eig$values, 1e-8)) %*% t(eig$vectors)
      cov_mat <- (cov_mat + t(cov_mat)) / 2
    }
    cov_mat
  }

  # Initialize Parameters (Starting Values)
  sigmaSq <- 1 / rgamma(1, shape = 1.5, rate = 0.5)
  a       <- switch(
    hyper_params$prior_type_a,
    gamma = rgamma(1, shape = hyper_params$s_a, rate = hyper_params$r_a),
    hcauchy = abs(rcauchy(1, location = 0, scale = hyper_params$scale_a))
  )
  b       <- switch(
    hyper_params$prior_type_b,
    gamma = rgamma(1, shape = hyper_params$s_b, rate = hyper_params$r_b),
    hcauchy = abs(rcauchy(1, location = 0, scale = hyper_params$scale_b))
  )
  phi     <- hyper_params$scale_phi
  if (hyper_params$prior_type_phi == "sqrt_hcauchy") phi <- phi^2
  log_nu <- log(rgamma(p, shape = 1))
  log_lambda <- log(rgamma(p, shape = 1))
  log_var <- log(phi) + log(a) - log(b) + log_nu - log_lambda
  beta <- rnorm(p, sd = exp(log_var/2))
  log_abs_beta <- log(abs(beta))
  if (sampler == "conditional") log_xi <- tpb_log_gamma(p, a) + log_nu

  # Storage matrices
  store_beta    <- matrix(0, nrow = n_save, ncol = p)
  store_beta_loglik <- store_total_logpost <-
    if (sampler == "marginal") numeric(n_save) else NULL
  store_scalars <- matrix(0, nrow = n_save, ncol = 4) # sigmaSq, phi, a, b
  colnames(store_scalars) <- c("sigmaSq", "phi", "a", "b")
  idx <- 0
  beta_loglik <- NA_real_
  total_logpost <- NA_real_

  cat(sprintf("Starting %s Gibbs Sampler...\n", sampler))

  accept_count_a_phi <- 0
  accept_count_b <- 0

  if (sampler == "marginal") {
    scale_factor <- (2.4^2) / 2
    emp_cov <- diag(2) * 0.01
    current_proposal_cov <- regularize_log_proposal_cov(scale_factor * emp_cov)

    block_samples <- matrix(0, nrow = adapt_block_size, ncol = 2)
    block_accepts <- 0
    block_idx <- 1
  }

  for (iter in 1:total_iter) {
    save_draw <- iter > num_burnin && (iter - num_burnin) %% thin == 0
    gaussian <- Gibbs_beta(X = X, y = y, sigmaSq = sigmaSq,
                           log_local_var = log_var, return_log = TRUE,
                           woodbury = woodbury, diagX = diagX)
    beta <- gaussian$beta
    log_abs_beta <- gaussian$log_abs_beta
    sigmaSq <- Gibbs_sigmaSq(n, X, y, beta)

    log_params_before <- NULL
    if (sampler == "conditional") {
      local <- tpb_conditional_locals(log_abs_beta, log_lambda, log_xi, a, b, phi)
      log_nu <- local$nu; log_lambda <- local$lambda; log_xi <- local$xi
      target <- tpb_conditional_target(log_abs_beta, log_nu, log_lambda, log_xi, hyper_params)
      log_params <- log_params_before <- log(c(a, b, phi))
      for (j in 1:3) {
        density <- function(value) {
          proposal <- log_params
          proposal[j] <- value
          target(proposal)
        }
        log_params[j] <- tpb_slice(log_params[j], density, slice_width[j])
      }
      a <- exp(log_params[1]); b <- exp(log_params[2]); phi <- exp(log_params[3])
      log_var <- log_params[3] + log_params[1] - log_params[2] + log_nu - log_lambda
    } else {
      mh <- run_marginal_mh(beta, a, b, phi, hyper_params, current_proposal_cov,
                             max_log_mh_step, log_abs_beta)
      if (mh$accepted_a_phi) accept_count_a_phi <- accept_count_a_phi + 1
      if (mh$accepted_b) accept_count_b <- accept_count_b + 1
      a   <- mh$a
      b   <- mh$b
      phi <- mh$phi
      beta_loglik <- mh$log_lik
      total_logpost <- mh$total_logpost

      if (iter <= num_burnin) {
        block_samples[block_idx, ] <- c(log(a), log(phi))
        if (mh$accepted_a_phi) {block_accepts <- block_accepts + 1}

        if (block_idx == adapt_block_size) {
          current_r <- block_accepts / adapt_block_size
          log_scale <- log(scale_factor) + 0.1 * (current_r - r_opt)
          scale_factor <- exp(log_scale)
          emp_cov <- cov(block_samples) + diag(2) * 1e-6
          current_proposal_cov <- regularize_log_proposal_cov(scale_factor * emp_cov)

          block_idx <- 1
          block_accepts <- 0
        } else {
          block_idx <- block_idx + 1
        }
      }
    }

    # Marginal MH discards the old local scales. Restore their conditional
    # distribution under the updated a, b, phi before the next beta update.
    if (sampler == "marginal")
      log_var <- tpb_draw_local_var(beta,a,b,phi,log_abs_beta,return_log=TRUE)
    if (interweave && iter %% interweave_every == 0) {
      log_params_before <- log(c(a, b, phi))
      moved <- tpb_interweave(X,y,log_var,a,b,phi,sigmaSq,hyper_params,
                             woodbury,diagX,slice_width)
      a <- moved$a; b <- moved$b; phi <- moved$phi
      beta <- moved$beta; log_var <- moved$log_var
      log_abs_beta <- moved$log_abs_beta
      if (save_draw) {
        density <- log_marginal_posterior(a, b, phi, beta, hyper_params, log_abs_beta)
        beta_loglik <- density$log_lik
        total_logpost <- density$log_posterior
      }
    }

    if (adapt_slice && iter <= num_burnin && !is.null(log_params_before)) {
      slice_updates <- slice_updates + 1L
      slice_jumps <- slice_jumps + abs(log(c(a, b, phi)) - log_params_before)
      if (slice_updates %% adapt_block_size == 0) {
        width <- 2 * slice_jumps / adapt_block_size
        slice_width[width > 0] <- width[width > 0]
        slice_jumps[] <- 0
      }
    }
    if (save_draw) {
      idx <- idx + 1
      store_beta[idx, ] <- beta
      if (sampler == "marginal") {
        store_beta_loglik[idx] <- beta_loglik
        store_total_logpost[idx] <- total_logpost
      }
      store_scalars[idx, ] <- c(sigmaSq, phi, a, b)
    }
  }

  accept_a_phi <- if (sampler == "marginal") accept_count_a_phi / total_iter else NA_real_
  accept_b <- if (sampler == "marginal") accept_count_b / total_iter else NA_real_

  result <- list(
    samples = list(
      beta = store_beta,
      beta_loglik = store_beta_loglik,
      total_logpost = store_total_logpost,
      scalars = store_scalars
    ),
    acceptance_rates = list(
      a = accept_a_phi,
      b = accept_b,
      phi = accept_a_phi
    ),
    hyper_params = hyper_params,
    sampler = sampler,
    interweave_every = if (interweave) interweave_every else NULL,
    slice_width = setNames(slice_width, c("a", "b", "phi"))
  )
  if (sampler == "marginal") {
    result$final_proposal_cov <- current_proposal_cov
    result$final_scale_factor <- scale_factor
  }

  return(result)
}
