#' Fully Gibbs Sampler Function with Miller's Approximation
#'
#' @import mvtnorm
#' @import coda
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param num_output Integer, number of posterior samples to save (after burn-in)
#' @param num_burnin Integer, number of burn-in iterations to discard
#' @param woodbury Logical, use Woodbury identity in beta update
#' @param diagX Logical, assume diagonal X
#' @param proposal_type String, "separate" or "bi_fixed" or "bi_adaptive" or "all_adaptive"
#' @param mh_step_a, mh_step_b, mh_step_phi step size for each hyper-parameter under the 'separate' mode
#' @param adapt_block_size Integer, number of iterations per adaptation block
#' @param r_opt Numeric, target acceptance rate for adaptive MH (default 0.3)
#' @param cov_matrix Matrix, 2x2 covariance matrix for bi_fixed proposals
#' @param hyper_params List of hyperparameters (prior_type_a, prior_type_b, s_a, r_a, s_b, r_b, scale_a, scale_b, scale_phi)
#' @param mh_step scalar, step size for MH uniform proposal
#' @return A list containing posterior samples matrices, acceptance rates, covariance matrix.
#' @export
fullGibbs <- function(X, y, num_output = 10000, num_burnin = 10000, thin = 1,
                      woodbury = TRUE, diagX = FALSE, proposal_type = "separate",
                      mh_step_a = 0.1, mh_step_b = 0.1, mh_step_phi = 0.1,
                      adapt_block_size = 100, r_opt = 0.3,
                      cov_matrix = matrix(c(0.01, -0.005, -0.005, 0.01), 2, 2),
                      hyper_params = list(prior_type_a = "gamma", prior_type_b = "gamma",
                                          s_a = 1.5, r_a = 1, s_b = 1.5, r_b = 1,
                                          scale_a = 1, scale_b = 1, scale_phi = 1)) {

  n <- nrow(X)
  p <- ncol(X)
  total_iter <- num_output + num_burnin
  n_save <- floor(num_output / thin)
  hyper_params <- modifyList(
    list(prior_type_a = "gamma", prior_type_b = "gamma",
         s_a = 1.5, r_a = 1, s_b = 1.5, r_b = 1,
         scale_a = 1, scale_b = 1, scale_phi = 1),
    hyper_params
  )
  hyper_params$prior_type_a <- match.arg(hyper_params$prior_type_a, c("gamma", "hcauchy"))
  hyper_params$prior_type_b <- match.arg(hyper_params$prior_type_b, c("gamma", "hcauchy"))

  # Initialize Parameters (Starting Values)
  sigmaSq <- 1 / rgamma(1, shape = 1.5, rate = 0.5) # var(y)
  a       <- if (hyper_params$prior_type_a == "gamma") {
    rgamma(1, shape = hyper_params$s_a, rate = hyper_params$r_a)
  } else {
    abs(rcauchy(1, location = 0, scale = hyper_params$scale_a))
  }
  b       <- if (hyper_params$prior_type_b == "gamma") {
    rgamma(1, shape = hyper_params$s_b, rate = hyper_params$r_b)
  } else {
    abs(rcauchy(1, location = 0, scale = hyper_params$scale_b))
  }
  a       <- pmin(pmax(a, 1e-4), 1e4)
  b       <- pmin(pmax(b, 1e-4), 1e4)
  phi     <- abs(rcauchy(1, location = 0, scale = hyper_params$scale_phi))
  phi     <- pmin(pmax(phi, 1e-4), 1e4) # 1
  nu      <- rgamma(p, shape = a, rate = 1) # 1
  lambda  <- rgamma(p, shape = b, rate = 1) # 1
  beta    <- rnorm(p, mean = 0, sd = sqrt(phi * (nu / lambda)))

  # w       <- 1
  # xi      <- rep(1, p)

  # Storage matrices (only for saved samples)
  store_beta    <- matrix(0, nrow = n_save, ncol = p)
  store_beta_loglik <- numeric(n_save)
  store_total_logpost <- numeric(n_save)
  store_scalars <- matrix(0, nrow = n_save, ncol = 4) # sigmaSq, phi, a, b; if use w then set it to 5
  colnames(store_scalars) <- c("sigmaSq", "phi", "a", "b")
  idx <- 0
  beta_loglik <- NA_real_
  total_logpost <- NA_real_

  # Progress bar
  #pb <- txtProgressBar(min = 0, max = total_iter, style = 3)

  cat(sprintf("Starting Marginalized Gibbs Sampler (Mode: %s)...\n", proposal_type))

  accept_count_a <- 0
  accept_count_b <- 0
  accept_count_phi <- 0

  if (proposal_type == "bi_adaptive") {
    d <- 2
    scale_factor <- (2.4^2) / d
    # emp_cov <- cov_matrix / scale_factor
    emp_cov <- diag(d) * 0.01
    current_proposal_cov <- scale_factor * emp_cov

    block_samples <- matrix(0, nrow = adapt_block_size, ncol = d)
    block_accepts <- 0
    block_idx <- 1
  }

  if (proposal_type == "all_adaptive") {
    d <- 3
    scale_factor <- (2.4^2) / d
    emp_cov <- diag(d) * 0.01
    current_proposal_cov <- scale_factor * emp_cov

    block_samples <- matrix(0, nrow = adapt_block_size, ncol = d)
    block_accepts <- 0
    block_idx <- 1
  }


  for (iter in 1:total_iter) {
    beta <- Gibbs_beta(X = X, y = y, a = a, b = b,
                       phi = phi, sigmaSq = sigmaSq, nu = nu, lambda = lambda,
                       woodbury = woodbury, diagX = diagX)

    sigmaSq <- Gibbs_sigmaSq(n, X, y, beta)
    lambda <- Gibbs_lambda(b, p, phi, beta, nu)
    #w <- Gibbs_w(phi)
    #phi <- Gibbs_phi(p, nu, lambda, w, beta)
    #xi <- Gibbs_xi(a, p, nu)
    nu <- Gibbs_nu(phi, beta, lambda, a)

    if (proposal_type == "separate"){
      a_new <- run_marginal_mh_uni(beta_vec = beta, target_param = "a",
                                   current_a = a, current_b = b, current_phi = phi,
                                   s_prior_a = hyper_params$s_a, r_prior_a = hyper_params$r_a,
                                   s_prior_b = hyper_params$s_b, r_prior_b = hyper_params$r_b,
                                   scale_phi = hyper_params$scale_phi,
                                   prior_type_a = hyper_params$prior_type_a,
                                   prior_type_b = hyper_params$prior_type_b,
                                   scale_a = hyper_params$scale_a,
                                   scale_b = hyper_params$scale_b,
                                   mh_step = mh_step_a)
      if (a_new$accepted) {accept_count_a <- accept_count_a + 1}
      a <- a_new$value
      beta_loglik <- a_new$log_lik
      total_logpost <- a_new$total_logpost

      b_new <- run_marginal_mh_uni(beta_vec = beta, target_param = "b",
                                   current_a = a, current_b = b, current_phi = phi,
                                   s_prior_a = hyper_params$s_a, r_prior_a = hyper_params$r_a,
                                   s_prior_b = hyper_params$s_b, r_prior_b = hyper_params$r_b,
                                   scale_phi = hyper_params$scale_phi,
                                   prior_type_a = hyper_params$prior_type_a,
                                   prior_type_b = hyper_params$prior_type_b,
                                   scale_a = hyper_params$scale_a,
                                   scale_b = hyper_params$scale_b,
                                   mh_step = mh_step_b)
      if (b_new$accepted) {accept_count_b <- accept_count_b + 1}
      b <- b_new$value
      beta_loglik <- b_new$log_lik
      total_logpost <- b_new$total_logpost

      phi_new <- run_marginal_mh_uni(beta_vec = beta, target_param = "phi",
                                     current_a = a, current_b = b, current_phi = phi,
                                     s_prior_a = hyper_params$s_a, r_prior_a = hyper_params$r_a,
                                     s_prior_b = hyper_params$s_b, r_prior_b = hyper_params$r_b,
                                     scale_phi = hyper_params$scale_phi,
                                     prior_type_a = hyper_params$prior_type_a,
                                     prior_type_b = hyper_params$prior_type_b,
                                     scale_a = hyper_params$scale_a,
                                     scale_b = hyper_params$scale_b,
                                     mh_step = mh_step_phi)
      if (phi_new$accepted) {accept_count_phi <- accept_count_phi + 1}
      phi <- phi_new$value
      beta_loglik <- phi_new$log_lik
      total_logpost <- phi_new$total_logpost
    } else if (proposal_type %in% c("bi_fixed", "bi_adaptive")) {
      prop_cov <- if(proposal_type == "bi_fixed") cov_matrix else current_proposal_cov
      a_phi_new <- run_marginal_mh_bi_a_phi(beta_vec = beta, current_a = a, current_b = b, current_phi = phi,
                                            s_prior_a = hyper_params$s_a, r_prior_a = hyper_params$r_a,
                                            s_prior_b = hyper_params$s_b, r_prior_b = hyper_params$r_b,
                                            scale_phi = hyper_params$scale_phi,
                                            prior_type_a = hyper_params$prior_type_a,
                                            prior_type_b = hyper_params$prior_type_b,
                                            scale_a = hyper_params$scale_a,
                                            scale_b = hyper_params$scale_b,
                                            cov_matrix = prop_cov)
      if (a_phi_new$accepted) {
        accept_count_a <- accept_count_a + 1
        accept_count_phi <- accept_count_phi + 1
      }
      a   <- a_phi_new$a
      phi <- a_phi_new$phi
      beta_loglik <- a_phi_new$log_lik
      total_logpost <- a_phi_new$total_logpost

      if (proposal_type == "bi_adaptive" && iter <= num_burnin) {
        block_samples[block_idx, ] <- c(log(a), log(phi))
        if (a_phi_new$accepted) {block_accepts <- block_accepts + 1}

        if (block_idx == adapt_block_size) {
          current_r <- block_accepts / adapt_block_size
          # learning rate = 0.1
          log_scale <- log(scale_factor) + 0.1 * (current_r - r_opt)
          scale_factor <- exp(log_scale)
          emp_cov <- cov(block_samples) + diag(2) * 1e-6
          current_proposal_cov <- scale_factor * emp_cov

          block_idx <- 1
          block_accepts <- 0
        } else {
          block_idx <- block_idx + 1
        }
      }

      b_new <- run_marginal_mh_uni(beta_vec = beta, target_param = "b",
                                   current_a = a, current_b = b, current_phi = phi,
                                   s_prior_a = hyper_params$s_a, r_prior_a = hyper_params$r_a,
                                   s_prior_b = hyper_params$s_b, r_prior_b = hyper_params$r_b,
                                   scale_phi = hyper_params$scale_phi,
                                   prior_type_a = hyper_params$prior_type_a,
                                   prior_type_b = hyper_params$prior_type_b,
                                   scale_a = hyper_params$scale_a,
                                   scale_b = hyper_params$scale_b,
                                   mh_step = mh_step_b)
      if (b_new$accepted) {accept_count_b <- accept_count_b + 1}
      b <- b_new$value
      beta_loglik <- b_new$log_lik
      total_logpost <- b_new$total_logpost
    } else if (proposal_type == "all_adaptive") {
      a_b_phi_new <- run_marginal_mh_tri_a_b_phi(beta_vec = beta, current_a = a, current_b = b, current_phi = phi,
                                                 s_prior_a = hyper_params$s_a, r_prior_a = hyper_params$r_a,
                                                 s_prior_b = hyper_params$s_b, r_prior_b = hyper_params$r_b,
                                                 scale_phi = hyper_params$scale_phi,
                                                 prior_type_a = hyper_params$prior_type_a,
                                                 prior_type_b = hyper_params$prior_type_b,
                                                 scale_a = hyper_params$scale_a,
                                                 scale_b = hyper_params$scale_b,
                                                 cov_matrix = current_proposal_cov)
      if (a_b_phi_new$accepted) {
        accept_count_a <- accept_count_a + 1
        accept_count_b <- accept_count_b + 1
        accept_count_phi <- accept_count_phi + 1
      }
      a   <- a_b_phi_new$a
      b   <- a_b_phi_new$b
      phi <- a_b_phi_new$phi
      beta_loglik <- a_b_phi_new$log_lik
      total_logpost <- a_b_phi_new$total_logpost

      if (iter <= num_burnin) {
        block_samples[block_idx, ] <- c(log(a), log(b), log(phi))
        if (a_b_phi_new$accepted) {block_accepts <- block_accepts + 1}

        if (block_idx == adapt_block_size) {
          current_r <- block_accepts / adapt_block_size
          log_scale <- log(scale_factor) + 0.1 * (current_r - r_opt)
          scale_factor <- exp(log_scale)
          emp_cov <- cov(block_samples) + diag(3) * 1e-6
          current_proposal_cov <- scale_factor * emp_cov

          block_idx <- 1
          block_accepts <- 0
        } else {
          block_idx <- block_idx + 1
        }
      }
    } else {
      stop("Invalid proposal_type! Choose 'separate', 'bi_fixed', 'bi_adaptive', or 'all_adaptive'.")
    }

    if (iter > num_burnin) {
      if ((iter - num_burnin) %% thin == 0) {
        idx <- idx + 1
        if (idx <= n_save) {
          store_beta[idx, ] <- beta
          store_beta_loglik[idx] <- beta_loglik
          store_total_logpost[idx] <- total_logpost
          store_scalars[idx, ] <- c(sigmaSq, phi, a, b)
        }
      }
    }
    # Update progress
    # if (iter %% 100 == 0) setTxtProgressBar(pb, iter)
  }

  accept_a <- accept_count_a / total_iter
  accept_b <- accept_count_b / total_iter
  accept_phi <- accept_count_phi / total_iter
  # close(pb)

  ### Diagnostics & Summary

  scalar_names <- colnames(store_scalars)
  result <- list(
    samples = list(
      beta = store_beta,
      beta_loglik = store_beta_loglik,
      total_logpost = store_total_logpost,
      scalars = store_scalars
    ),
    acceptance_rates = list(
      a = accept_a,
      b = accept_b,
      phi = accept_phi
    ),
    hyper_params = hyper_params
  )
  if (proposal_type %in% c("bi_adaptive", "all_adaptive")) {
    result$final_proposal_cov <- current_proposal_cov
    result$final_scale_factor <- scale_factor
  }

  return(result)
}
