#' @param beta_samples S*p matrix, posterior beta samples
#' @param sigmaSq numeric, estimated sigmaSq
#' @return log-likelihood matrix, (i, j)-entry is the log-likelihood of y_i under the j-th sampled beta.
log_lik <- function(X, y, beta_samples, sigmaSq) {
  pred_mat <- X %*% t(beta_samples)
  y_mat <- matrix(rep(y, ncol(pred_mat)), nrow = nrow(X)) # y_mat: n*S matrix, each column is a replication of y
  log_lik_mat <- dnorm(y_mat, mean = pred_mat, sd = sqrt(sigmaSq), log = TRUE)
  return(log_lik_mat)
}

#' Wrapper function to calculate the score for a single initialization.
#' @param beta_k posterior mean after a given initialization and moderate number of iterations.
#' @param sigmaSq_k estimated sigmaSq after a given initialization and moderate number of iterations.
#' @return numeric, observed-data log-likelihood
get_initialization_score <- function(X, y, beta_k, sigmaSq_k) {
  # Converts a point estimate beta_k into a 1-row matrix to be used
  # by the user's log_lik function, then sums the result to get a single score.
  
  beta_k_mtx <- matrix(beta_k, nrow = 1)
  log_lik_matrix <- log_lik(X = X, y = y, beta_samples = beta_k_mtx, sigmaSq = sigmaSq_k)
  log_lik_score <- sum(log_lik_matrix)
  return(log_lik_score)
}


#' Calculate the Complete-Data Log-Likelihood for a Given Model
#'
#' This function computes the logarithm of p(y|beta, model_k)xp(beta|lambda, nu, model_k)xp(lambda|model_k)xp(nu|model_k)
#' which serves as the weight in the model selection step.
#' @param beta_vec, nu_vec, lambda_vec The current p-dimensional vectors for the latent variables.
#' @param X, y The data.
#' @param model_params A list containing the specific hyperparameters (a, b, omega, sigmaSq) for a candidate model.
#' @return A single numeric value representing the complete-data log-likelihood.
calculate_complete_loglik <- function(beta_vec, nu_vec, lambda_vec, X, y, model_params) {
  a <- model_params$a
  b <- model_params$b
  omega <- model_params$omega
  sigmaSq <- model_params$sigmaSq
  
  log_lik_y <- sum(dnorm(y, mean = X %*% beta_vec, sd = sqrt(sigmaSq), log = TRUE))
  log_prior_beta <- sum(dnorm(beta_vec, mean = 0, sd = sqrt(omega * nu_vec / (lambda_vec + .Machine$double.eps)), log = TRUE))
  log_prior_nu <- sum(dgamma(nu_vec + .Machine$double.eps, shape = a, rate = a, log = TRUE))
  log_prior_lambda <- sum(dgamma(lambda_vec + .Machine$double.eps, shape = b, rate = b, log = TRUE))
  
  total_log_lik <- log_lik_y + log_prior_beta + log_prior_nu + log_prior_lambda
  
  return(c(total = total_log_lik, 
           y = log_lik_y,
           beta = log_prior_beta,
           nu = log_prior_nu,
           lambda = log_prior_lambda
  ))
}


#' @return
calculate_marginal_loglik_beta <- function(beta_vec, model_params) {
  
  d_tpb <- function(beta_vec, a, b, phi) {
    
    tricomi_U <- function(a, b, z) {
      integrand <- function(t) {
        exp(-z * t) * t^(a - 1) * (1 + t)^(b - a - 1)
      }
      integral_result <- tryCatch({
        integrate(integrand, lower = 0, upper = Inf)$value
      }, error = function(e) { return(NA) })
      return(integral_result / gamma(a))
    }
    
    sapply(beta_vec, function(beta_i) {
      const <- gamma(0.5 + b) * gamma(a + b) / (gamma(a) * gamma(b) * sqrt(2 * pi * phi))
      z_val <- beta_i^2 / (2 * phi)
      U_val <- tricomi_U(0.5 + b, 1.5 - a, z_val)
      if (is.na(U_val)) return(NA)
      return(const * U_val)
    })
  }
  
  
  a <- model_params$a
  b <- model_params$b
  omega <- model_params$omega
  phi <- omega * b / a
  loglik_individual <- log(pmax(d_tpb(beta_vec = beta_vec, a = a, b = b, phi = phi), .Machine$double.xmin))
  
  if (any(!is.finite(loglik_individual))) {
    return(-Inf)
  }
  
  return(sum(loglik_individual))
}

