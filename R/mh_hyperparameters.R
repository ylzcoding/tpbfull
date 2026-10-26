# Joint MH for log(a), log(phi), followed by a scalar update of log(b).
run_marginal_mh <- function(beta_vec, a, b, phi, hyper_params, cov_matrix,
                            max_log_mh_step = Inf,
                            log_abs_beta = log(abs(beta_vec))) {
  log_density <- function(param) {
    log_marginal_posterior(param[1], param[2], param[3], beta_vec,
                           hyper_params, log_abs_beta)
  }
  current <- c(a, b, phi)
  log_cur <- log(current[c(1, 3)])
  log_new <- as.vector(mvtnorm::rmvnorm(1, mean = log_cur, sigma = cov_matrix))
  current_density <- log_density(current)
  accepted_a_phi <- FALSE
  if (!(is.finite(max_log_mh_step) && any(abs(log_new-log_cur) > max_log_mh_step))) {
    proposal <- current
    proposal[c(1, 3)] <- exp(log_new)
    proposal_density <- log_density(proposal)
    log_r <- proposal_density$log_posterior - current_density$log_posterior +
      (log_new[1]-log_cur[1]) + (log_new[2]-log_cur[2])
    if (!is.na(log_r) && log(runif(1)) < log_r) {
      current <- proposal
      current_density <- proposal_density
      accepted_a_phi <- TRUE
    }
  }

  log_b <- log(current[2])
  log_b_new <- rnorm(1, mean = log_b, sd = 0.1)
  accepted_b <- FALSE
  if (!(is.finite(max_log_mh_step) && abs(log_b_new-log_b) > max_log_mh_step)) {
    proposal <- current
    proposal[2] <- exp(log_b_new)
    proposal_density <- log_density(proposal)
    log_r <- proposal_density$log_posterior - current_density$log_posterior +
      log(proposal[2]) - log(current[2])
    if (!is.na(log_r) && log(runif(1)) < log_r) {
      current <- proposal
      current_density <- proposal_density
      accepted_b <- TRUE
    }
  }
  list(a = current[1], b = current[2], phi = current[3],
       accepted_a_phi = accepted_a_phi, accepted_b = accepted_b,
       log_lik = current_density$log_lik, total_logpost = current_density$log_posterior)
}
