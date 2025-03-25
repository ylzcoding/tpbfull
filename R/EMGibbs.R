#' @import mvtnorm
#' @param iter_burnin number of burn-in samples in each iteration
#' @param num_output number of samples generated after convergence
#' @param num_burnin number of samples burnin after convergence
#' @param max_iter numeric, max number of iterations, defaults to 0 (just get samples at starting values or provided values)
#' @param fix.a numeric, fix the value of a to skip estimation
#' @param omega_hyper numeric, hyper-parameter of the hyper-prior of omega, omega^{1/2} ~ half-cauchy(0, omega_hyper^{1/2})
#' @param a.start starting point of a for EM-within-Gibbs
#' @param delta1 numeric, a small number added to avoid division by zero
#' @param delta2 numeric, error bound for convergence in isConverged_hyper
#' @param window_size numeric, size of the sliding window to determine the convergence
#' @param option "augmented" or "naive", specifying how we would like to sample nu
#' @return list of posterior samples
#' @export
EM_GibbsTPB = function(num_output, num_burnin, iter_burnin, X, y, max_iter = 0,
                       fix.a = NULL, fix.b = NULL, fix.omega = NULL,
                       omega_hyper = NULL, a.start = 0.5, b.start = 0.5, omega.start = 1,
                       delta1 = 1e-3, delta2 = 1e-3, window_size = 5, option = "augmented") {
  # if fix.a is not NULL, fix a at that value
  n = dim(X)[1]
  p = dim(X)[2]
  emp_omega = is.null(omega_hyper) # is omega_hyper = NULL, then estimate omega using EM-within-Gibbs or fix omega

  a = ifelse(is.null(fix.a), a.start, fix.a)
  b = ifelse(is.null(fix.b), b.start, fix.b)

  a_vec = c(a)
  b_vec = c(b)

  if (emp_omega) {
    omega = ifelse(is.null(fix.omega), omega.start, fix.omega)
    omega_vec = c(omega)
  } else{
    omega = rgamma(1, 0.5, 0.1)
    w = rgamma(1, 0.5, 0.01)
  }

  ###### Initialization
  sigmaSq0 = var(y)
  beta0 = t(mvtnorm::rmvnorm(1, sigma = sigmaSq0 * diag(p)))
  nu0 = rgamma(p, a, 1)
  lambda0 = rgamma(p, b, 1)
  if (option == "augmented") {
    xi0 = a*sqrt(exp(log(beta0^2) + log(lambda0) - log(sigmaSq0) - log(omega))/2)
  } else {
    xi0 = NULL
  }
  ######
  if (max_iter > 0) {
    for (k in 1:max_iter) {
      iter_sample = ifelse(k < 40, 1000, ifelse(k < 100, 5000, 10000))
      if (emp_omega) {
        samples = getsamples.emp_omega(num_output = iter_sample, num_burnin = iter_burnin,
                                       X = X, y = y, a = a_vec[k], b = b_vec[k], omega = omega_vec[k],
                                       sigmaSq0 = sigmaSq0, beta0 = beta0, nu0 = nu0, lambda0 = lambda0, xi0 = xi0)
      } else {
        samples = getsamples.full_omega(num_output = iter_sample, num_burnin = iter_burnin,
                                        X = X, y = y, a = a_vec[k], b = b_vec[k], omega_hyper = omega_hyper,
                                        omega0 = omega, w0 = w, sigmaSq0 = sigmaSq0,
                                        beta0 = beta0, nu0 = nu0, lambda0 = lambda0, xi0 = xi0)
      }

      ###### M-step
      if (is.null(fix.a)) {
        a_vec = c(a_vec, M.step(samples$nu))
      } else {
        a_vec = c(a_vec, a) # if fix.a is not NULL, skip estimating a
      }

      if (is.null(fix.b)) {
        b_vec = c(b_vec, M.step(samples$lambda))
      } else {
        b_vec = c(b_vec, b) # if fix.b is not NULL, skip estimating a
      }

      if (emp_omega) {
        if (is.null(fix.omega)) {
          omega_vec = c(omega_vec, M.step_omega(samples$beta, samples$lambda, samples$nu, samples$sigmaSq))
        } else { # if fix.omega is not NULL, skip estimating omega
          omega_vec = c(omega_vec, omega)
        }
      }

      ###### reset initial values for the next iteration
      beta0 = samples$beta[iter_sample, ]
      nu0 = samples$nu[iter_sample, ]
      lambda0 = samples$lambda[iter_sample, ]
      xi0 = samples$xi[iter_sample, ]
      sigmaSq0 = samples$sigmaSq[iter_sample]
      if (!emp_omega) {
        # if go with fully bayes on omega, we also need to reset initial values of omega and w
        omega = samples$omega[iter_sample]
        w = samples$w[iter_sample]
      }

      ###### convergence?
      if (emp_omega) {
        if (isConverged_hyper(a_vec, b_vec, omega_vec, delta1, delta2, window_size)) {
          cat("Convergence achieved after", k, "iterations.\n")
          break
        }
        cat(k, "-th iteration completed.\n")
        cat("a=", a_vec[k+1], "\n")
        cat("b=", b_vec[k+1], "\n")
        cat("omega=", ifelse(emp_omega, omega_vec[k+1], mean(samples$omega)), "\n")
      }
    }
  }

  ### get posterior samples after convergence
  l = length(a_vec)
  a.est = a_vec[l]
  b.est = b_vec[l]
  if (emp_omega) {
    omega.est = omega_vec[l]
    output = getsamples.emp_omega(num_output = num_output, num_burnin = num_burnin,
                                  X = X, y = y, a = a.est, b = b.est, omega = omega.est,
                                  sigmaSq0 = sigmaSq0, beta0 = beta0, nu0 = nu0, lambda0 = lambda0, xi0 = xi0)
  } else {
    output = getsamples.full_omega(num_output = num_output, num_burnin = num_burnin,
                                   X = X, y = y, a = a.est, b = b.est, omega_hyper = omega_hyper,
                                   omega0 = omega, w0 = w, sigmaSq0 = sigmaSq0,
                                   beta0 = beta0, nu0 = nu0, lambda0 = lambda0, xi0 = xi0)
  }

  res = list(a = a_vec, b = b_vec, beta = output$beta, sigmaSq = output$sigmaSq)
  if (emp_omega) {
    res$omega <- omega_vec
  } else {
    res$omega <- output$omega
  }
  return(res)
}
