#' @return logarithm of the joint prior, likelihood, and joint posterior
#' @export
log_res <- function(a, b, omega, sigmaSq, y, X, beta, lambda, nu) {
  log_prior <-
    sum(
      dgamma(lambda, shape = a, rate = a, log = TRUE)
      + dgamma(nu, shape = b, rate = b, log = TRUE)
      + dnorm(beta, mean = 0, sd = sqrt(omega * nu / lambda), log = TRUE)
    )
  log_lik <- sum(dnorm(y, mean = X %*% beta, sd = sqrt(sigmaSq), log = TRUE))
  log_post <- log_prior + log_lik

  return(
    list(log_post = log_post,
         log_prior = log_prior,
         log_lik = log_lik)
  )
}

#' @param beta_mat M by p sample matrix
#' @return Important sampling weights (w_1^k, ..., w_M^k)
#' @export
IS_weights <- function(y, X, beta_mat, lambda_mat, nu_mat,
                       a_cur, b_cur, omega_cur, sigmaSq_cur,
                       a0, b0, omega0, sigmaSq0) {
  ## compute w_i^k = L(y, beta(i)^0, lambda(i)^0, lambda(i)^0 | a^k, b^k, omega^k, sigmaSq^k) /
  ##                 L(y, beta(i)^0, lambda(i)^0, lambda(i)^0 | a^0, b^0, omega^0, sigmaSq^0)
  ## where beta(i)^0 is the i-th beta vector sampled using a^0, b^0, omega^0, sigmaSq^0
  M <- nrow(beta_mat)
  log_post_cur <- vapply(1:M, function(i) {
    log_res(a = a_cur,
            b = b_cur,
            omega = omega_cur,
            sigmaSq = sigmaSq_cur,
            y = y,
            X = X,
            beta = beta_mat[i, ],
            lambda = lambda_mat[i, ],
            nu = nu_mat[i, ]
            )$log_post
  }, numeric(1))

  log_post0 <- vapply(1:M, function(i) {
    log_res(a = a0,
            b = b0,
            omega = omega0,
            sigmaSq = sigmaSq0,
            y = y,
            X = X,
            beta = beta_mat[i, ],
            lambda = lambda_mat[i, ],
            nu = nu_mat[i, ]
    )$log_post
  }, numeric(1))
  
  logw <- log_post_cur - log_post0
  logw <- logw - max(logw)
  IS_weights <- exp(logw)
  return(IS_weights)
}






