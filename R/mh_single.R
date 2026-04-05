#' single variable Marginal MH
#' @param target_param ("a", "b", or "phi")
run_marginal_mh_uni <- function(beta_vec, target_param,
                                current_a, current_b, current_phi,
                                s_prior_a, r_prior_a,
                                s_prior_b, r_prior_b,
                                scale_phi, mh_step = 0.1) {

  current_val <- switch(target_param,
                        "a" = current_a,
                        "b" = current_b,
                        "phi" = current_phi)

  new_val <- exp(rnorm(1, mean = log(current_val), sd = mh_step))

  a_new <- ifelse(target_param == "a", new_val, current_a)
  b_new <- ifelse(target_param == "b", new_val, current_b)
  phi_new <- ifelse(target_param == "phi", new_val, current_phi)

  log_post_old <- log_marginal_posterior(current_a, current_b, current_phi, beta_vec,
                                         s_prior_a, r_prior_a, s_prior_b, r_prior_b, scale_phi)
  log_post_new <- log_marginal_posterior(a_new, b_new, phi_new, beta_vec,
                                         s_prior_a, r_prior_a, s_prior_b, r_prior_b, scale_phi)

  log_r <- log_post_new - log_post_old + log(new_val) - log(current_val)

  if (!is.na(log_r) && log(runif(1)) < log_r) {
    return(new_val)
  } else {
    return(current_val)
  }
}



run_marginal_mh_bi_a_phi <- function(beta_vec, current_a, current_b, current_phi,
                                     s_prior_a, r_prior_a, s_prior_b, r_prior_b, scale_phi,
                                     cov_matrix = matrix(c(0.01, -0.005, -0.005, 0.01), 2, 2)) {

  log_cur <- c(log(current_a), log(current_phi))
  log_new <- as.vector(mvtnorm::rmvnorm(1, mean = log_cur, sigma = cov_matrix))

  a_new   <- exp(log_new[1])
  phi_new <- exp(log_new[2])

  log_post_old <- log_marginal_posterior(current_a, current_b, current_phi, beta_vec,
                                         s_prior_a, r_prior_a, s_prior_b, r_prior_b, scale_phi)
  log_post_new <- log_marginal_posterior(a_new, current_b, phi_new, beta_vec,
                                         s_prior_a, r_prior_a, s_prior_b, r_prior_b, scale_phi)

  log_r <- log_post_new - log_post_old + (log_new[1] - log_cur[1]) + (log_new[2] - log_cur[2])

  if (!is.na(log_r) && log(runif(1)) < log_r) {
    return(list(a = a_new, phi = phi_new, accepted = TRUE))
  } else {
    return(list(a = current_a, phi = current_phi, accepted = FALSE))
  }
}

run_marginal_mh_tri_a_b_phi <- function(beta_vec, current_a, current_b, current_phi,
                                        s_prior_a, r_prior_a, s_prior_b, r_prior_b, scale_phi,
                                        cov_matrix = diag(3) * 0.01) {

  log_cur <- c(log(current_a), log(current_b), log(current_phi))
  log_new <- as.vector(mvtnorm::rmvnorm(1, mean = log_cur, sigma = cov_matrix))

  a_new   <- exp(log_new[1])
  b_new   <- exp(log_new[2])
  phi_new <- exp(log_new[3])

  log_post_old <- log_marginal_posterior(current_a, current_b, current_phi, beta_vec,
                                         s_prior_a, r_prior_a, s_prior_b, r_prior_b, scale_phi)
  log_post_new <- log_marginal_posterior(a_new, b_new, phi_new, beta_vec,
                                         s_prior_a, r_prior_a, s_prior_b, r_prior_b, scale_phi)

  log_r <- log_post_new - log_post_old + (log_new[1] - log_cur[1]) + (log_new[2] - log_cur[2]) + (log_new[3] - log_cur[3])

  if (!is.na(log_r) && log(runif(1)) < log_r) {
    return(list(a = a_new, b = b_new, phi = phi_new, accepted = TRUE))
  } else {
    return(list(a = current_a, b = current_b, phi = current_phi, accepted = FALSE))
  }
}




