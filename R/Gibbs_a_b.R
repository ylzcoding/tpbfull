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

#' Sample posterior for hyperparameter 'a'
#'
#' @param nu p*1 vector of local shrinkage parameters
#' @param a_current scalar, current value of 'a'
#' @param s_a scalar, shape parameter for prior on 'a'
#' @param r_a scalar, rate parameter for prior on 'a'
#' @return scalar a_new
#' @export
Gibbs_a <- function(nu, a, s_a = 1, r_a = 1) {
  # Apply Miller's method
  # Sample from the approximated Gamma(A, B) posterior
  params <- run_miller_algorithm(x = nu, mu = a, s_prior = s_a, r_prior = r_a)
  a_new <- rgamma(1, shape = params$A, rate = params$B)
  return(a_new)
}

#' Sample posterior for hyperparameter 'b'
#'
#' @param lambda p*1 vector of local shrinkage parameters
#' @param b_current scalar, current value of 'b'
#' @param s_b_prior scalar, shape parameter for prior on 'b' (renamed to avoid conflict with local b)
#' @param r_b_prior scalar, rate parameter for prior on 'b'
#' @return scalar b_new
#' @export
Gibbs_b <- function(lambda, b, s_b = 1, r_b = 1) {
  # Apply Miller's method
  # Sample from the approximated Gamma(A, B) posterior
  params <- run_miller_algorithm(x = lambda, mu = b, s_prior = s_b, r_prior = r_b)
  b_new <- rgamma(1, shape = params$A, rate = params$B)
  return(b_new)
}
