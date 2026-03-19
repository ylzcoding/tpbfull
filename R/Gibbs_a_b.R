#' Miller's Approximation Algorithm for Gamma Shape Posterior
#' Find the parameters A and B of a Gamma(A, B) distribution that approximates the posterior of the shape parameter 'a'.
#' Model: X ~ Gamma(a, a/mu) (implying mean = mu)
#' Prior: a ~ Gamma(s_prior, r_prior)
#'
#' @param x Vector of data (e.g., nu for parameter a, or lambda for parameter b)
#' @param mu Current value of the parameter (a or b). In the rate=1 parameterization, mu = a (or b).
#' @param s_prior Prior shape parameter for a (or b)
#' @param r_prior Prior rate parameter for a (or b)
#' @param tol Tolerance for convergence
#' @param max_iter Maximum number of iterations
#' @return A list containing A and B
run_miller_algorithm <- function(x, mu, s_prior, r_prior, tol = 1e-4, max_iter = 20) {

  n <- length(x)
  R <- sum(log(x))
  S <- sum(x)
  T <- S/mu - R + n * log(mu) - n

  # Initialization
  A <- s_prior + n / 2
  B <- r_prior + T

  # Iterative Update
  for (j in 1:max_iter) {
    a_est <- A / B # estimated mean
    A <- s_prior - n * a_est + n * (a_est^2) * trigamma(a_est)
    B <- r_prior + (A - s_prior) / a_est - n * log(a_est) + n * digamma(a_est) + T

    # Check Convergence
    if (abs(a_est / (A / B) - 1) < tol) {
      return(list(A = A, B = B))
    }
  }

  # Return best estimate if max_iter reached
  return(list(A = A, B = B))
}


#' Metropolis-Hastings Algorithm for Gamma Prior with log-normal proposal
#'
#' @param current_val (a or b)
#' @param x (nu or lambda)
#' @param s_prior
#' @param r_prior
#' @param w step size
#' @return
run_mh_gamma <- function(current_val, x, s_prior, r_prior, w = 0.05) {

  p <- length(x)
  sum_log_x <- sum(log(x))
  # log(new_val) ~ normal(log(current_val), w^2)
  new_val <- exp(rnorm(1, mean = log(current_val), sd = w))

  # log target distribution: g(val) = -n*log(Gamma(val)) + val*sum(log(x)) + (s-1)*log(val) - r*val
  log_post_new <- -p * lgamma(new_val) + new_val * sum_log_x + (s_prior - 1) * log(new_val) - r_prior * new_val
  log_post_old <- -p  * lgamma(current_val) + current_val * sum_log_x + (s_prior - 1) * log(current_val) - r_prior * current_val
  log_r <- log_post_new - log_post_old + log(new_val) - log(current_val)

  if (log(runif(1)) < log_r) {
    return(new_val)
  } else {
    return(current_val)
  }
}


#' Metropolis-Hastings Algorithm for Half-Cauchy Prior with log-normal proposal
#'
#' @param current_val (a or b)
#' @param x (nu or lambda)
#' @param scale_prior Half-Cauchy(0, scale)
#' @param w step_size
#' @return
run_mh_hcauchy <- function(current_val, x, scale_prior, w = 0.05) {

  n <- length(x)
  sum_log_x <- sum(log(x))

  new_val <- exp(rnorm(1, mean = log(current_val), sd = w))

  # log(Prior) proportional to -log(1 + (val / scale_prior)^2)
  log_prior_new <- -log(1 + (new_val / scale_prior)^2)
  log_prior_old <- -log(1 + (current_val / scale_prior)^2)

  # log(Likelihood) = -n * lgamma(val) + val * sum(log(x))
  log_lik_new <- -n * lgamma(new_val) + new_val * sum_log_x
  log_lik_old <- -n * lgamma(current_val) + current_val * sum_log_x

  log_post_new <- log_lik_new + log_prior_new
  log_post_old <- log_lik_old + log_prior_old
  log_r <- log_post_new - log_post_old + log(new_val) - log(current_val)

  if (log(runif(1)) < log_r) {
    return(new_val)
  } else {
    return(current_val)
  }
}



#' Sample posterior for hyperparameter 'a'
#'
#' @param nu p*1 vector of local shrinkage parameters
#' @param a scalar, current value of 'a'
#' @param prior string, "gamma" or "hcauchy"
#' @param s_a scalar, shape parameter for Gamma prior
#' @param r_a scalar, rate parameter for Gamma prior
#' @param scale_a scalar, scale parameter for Half-Cauchy prior
#' @param method string, "miller" or "mh" (only applicable if prior_option="gamma")
#' @param mh_step scalar, step size for MH uniform proposal
#' @return scalar a_new
#' @export
Gibbs_a <- function(nu, a, prior = "gamma",
                    s_a = 1.5, r_a = 1, scale_a = 1,
                    method = "mh", mh_step = 0.05) {

  if (prior == "gamma") {
    if (method == "miller") {
      params <- run_miller_algorithm(x = nu, mu = a, s_prior = s_a, r_prior = r_a)
      a_new <- rgamma(1, shape = params$A, rate = params$B)
      return(a_new)
    } else if (method == "mh") {
      a_new <- run_mh_gamma(current_val = a, x = nu, s_prior = s_a, r_prior = r_a, w = mh_step)
      return(a_new)
    } else {
      stop("Invalid method for gamma prior. Please choose 'miller' or 'mh'.")
    }

  } else if (prior == "hcauchy") {
    if (method == "miller") {
      warning("Miller's approximation is not applicable to Half-Cauchy prior. Falling back to 'mh'.")
    }
    a_new <- run_mh_hcauchy(current_val = a, x = nu, scale_prior = scale_a, w = mh_step)
    return(a_new)

  } else {
    stop("Invalid prior_option. Please choose 'gamma' or 'hcauchy'.")
  }
}


#' Sample posterior for hyperparameter 'b'
#'
#' @param lambda p*1 vector of local shrinkage parameters
#' @param b scalar, current value of 'b'
#' @param prior string, "gamma" or "hcauchy"
#' @param s_b scalar, shape parameter for Gamma prior
#' @param r_b scalar, rate parameter for Gamma prior
#' @param scale_b scalar, scale parameter for Half-Cauchy prior
#' @param method string, "miller" or "mh" (only applicable if prior_option="gamma")
#' @param mh_step scalar, step size for MH uniform proposal
#' @return scalar b_new
#' @export
Gibbs_b <- function(lambda, b, prior = "gamma",
                    s_b = 1.5, r_b = 1, scale_b = 1,
                    method = "mh", mh_step = 0.05) {

  if (prior == "gamma") {
    if (method == "miller") {
      params <- run_miller_algorithm(x = lambda, mu = b, s_prior = s_b, r_prior = r_b)
      b_new <- rgamma(1, shape = params$A, rate = params$B)
      return(b_new)
    } else if (method == "mh") {
      b_new <- run_mh_gamma(current_val = b, x = lambda, s_prior = s_b, r_prior = r_b, w = mh_step)
      return(b_new)
    } else {
      stop("Invalid method for gamma prior. Please choose 'miller' or 'mh'.")
    }

  } else if (prior == "hcauchy") {
    if (method == "miller") {
      warning("Miller's approximation is not applicable to Half-Cauchy prior. Falling back to 'mh'.")
    }
    b_new <- run_mh_hcauchy(current_val = b, x = lambda, scale_prior = scale_b, w = mh_step)
    return(b_new)

  } else {
    stop("Invalid prior_option. Please choose 'gamma' or 'hcauchy'.")
  }
}
