#' @import mvtnorm
#' @param iter_burnin number of burn-in samples in each iteration
#' @param iter_samples number of samples used for estimation in each iteration
#' @param output_burnin number of burn-in samples needed for the draw after convergence
#' @param output_samples number of samples produced after convergence
#' @param delta1 numeric, a small number added to avoid division by zero
#' @param delta2 numeric, error bound for convergence in isConverged_hyper
#' @param delta3 numeric, error bound for convergence in isConverged_rho
#' @param window_size numeric, size of the sliding window
#' @param converge_rule "rho" or "hyper", specifying what to monitor to decide the convergence
#' @param max_iter numeric, max number of iterations, defaults to 0 (just get samples at starting values or provided values)
#' @param IS_period number of iterations between each resampling step in importance sampling. Set to 1 by default to disable importance sampling
#' @param option "augmented" or "naive", specifying how we would like to sample nu
#' @param a.start starting point of a for EM-within-Gibbs if a.start = b.start = phi.start = sigmaSq.start = NULL (simply a.start = NULL), then run full Bayes algorithm with Metropolis-Hastings
#' @return list of posterior samples
#' @export
EM_GibbsTPB = function(iter_burnin, iter_samples, output_burnin, output_samples,
                       X, y, max_iter = 0, IS_period = 1, option = "augmented", 
                       omega_hyper = NULL, woodbury = FALSE, approx = FALSE,
                       delta2 = 1e-3, delta3 = 1e-3, delta1 = 1e-4, 
                       window_size = 5, converge_rule = "hyper",
                       a = NULL, b = NULL, omega = NULL, sigmaSq = NULL,
                       a.start = 0.5, b.start = 0.5, omega.start = 1/n, sigmaSq.start = 1,
                       tunea = 0.1, tuneb = 0.2, uppera = 1, upperb = 1, diagX = FALSE){

  # EM_step: number of samples used for empirical expectation estimation(in EM), could be a multiplier, such that the number of samples increases linearly with each iteration k, e.g. 500*k

  n = dim(X)[1]
  p = dim(X)[2]

  null.a = is.null(a)
  null.b = is.null(b)
  null.omega = is.null(omega)
  null.sigmaSq = is.null(sigmaSq)
  null.astart = is.null(a.start)


  if (null.astart){
    ### Initialize everything
    sigmaSq0 = var(y)
    a0 = runif(1, 0, 1)
    b0 = runif(1, 0, 1)
    beta0 = rnorm(p, 0, sigma = sqrt(sigmaSq0))
    nu0 = rgamma(p, a0, 1)
    lambda0 = rgamma(p, b0, 1)
    omega0 = rgamma(1, 0.5, 0.1)
    w0 = rgamma(1, 0.5, 0.01)
    if (option == "augmented") {
      xi0 = a0*sqrt(exp(log(beta0^2) + log(lambda0) - log(omega0))/2)
    } else {
      xi0 = NULL
    }

    return(fullGibbs(burnin_num = output_burnin, output_num = output_samples, X = X, y = y,
                     a0 = a0, b0 = b0, beta0 = beta0, nu0 = nu0, lambda0 = lambda0,
                     omega0 = omega0, w0 = w0, sigmaSq0 = sigmaSq0, omega_hyper = omega_hyper,
                     xi0 = xi0, tunea = tunea, tuneb = tuneb, uppera = uppera, upperb = upperb, 
                     woodbury = woodbury, approx = approx, diagX = diagX))
    } else {

    if (null.a) {a = a.start}
    if (null.b) {b = b.start}
    if (null.omega) {omega = omega.start}
    if (null.sigmaSq) {sigmaSq = sigmaSq.start}

    a_vec = c(a) # create vectors to record their trajectories
    b_vec = c(b)
    omega_vec = c(omega)
    sigmaSq_vec = c(sigmaSq)

    ### Initialization
    beta0 = rnorm(p, mean = 0, sd = sqrt(sigmaSq))
    nu0 = rgamma(p, 1, 1)
    lambda0 = rgamma(p, 1, 1)

    if (option == "augmented") {
      xi0 = a*sqrt(exp(log(beta0^2) + log(lambda0) - log(omega))/2)
    } else {
      xi0 = NULL
    }
    
    ### generate samples for EM algorithm and burn-in process
    if (max_iter > 0) {
      for (k in 1:max_iter){
        # num_sample = ifelse(adapt, k, 1)*EM_step
        # num_sample = ifelse(k < 40, 1000, ifelse(k < 100, 2500, 5000))
        num_sample = iter_samples
        if ((k-1) %% IS_period == 0) {
          # if we need resampling
          samples = getsamples.emp(
            num = num_sample, X = X, y = y,
            a = a_vec[k], b = b_vec[k], omega = omega_vec[k], sigmaSq = sigmaSq_vec[k],
            beta0 = beta0, nu0 = nu0, lambda0 = lambda0, burn = iter_burnin, 
            xi0 = xi0, woodbury = woodbury, approx = approx, diagX = diagX
            )
          
          # Regular M-step
          if (null.a) {
            a_vec = c(a_vec, M.step(samples$nu*a_vec[k]))
            } else {
              a_vec = c(a_vec, a)
            }
          
          if (null.b) {
            b_vec = c(b_vec, M.step(samples$lambda*b_vec[k]))
            } else {
              b_vec = c(b_vec, b)
            }
          
          if (null.sigmaSq) {
            sigmaSq_vec = c(sigmaSq_vec, M.step_sigmaSq(beta_matrix = samples$beta, X = X, y = y, n = n, diagX = diagX))
            } else {
              sigmaSq_vec = c(sigmaSq_vec, sigmaSq)
            }
          
          if (null.omega) {
            omega = M.step_phi(beta_matrix = samples$beta, lambda_matrix = samples$lambda, nu_matrix = samples$nu)
            omega_vec = c(omega_vec, omega)
            } else {
              omega_vec = c(omega_vec, omega)
            }
          
        # reset initial values for the next resampling
          beta0 = samples$beta[num_sample, ]
          nu0 = samples$nu[num_sample, ]
          lambda0 = samples$lambda[num_sample, ]
        # record this iteration number
          resample_iter = k
        } else {
          # if we use existed samples
          # M-step with importance sampling
          weights_vec = IS_weights(y = y, X = X,
                                   beta_mat = samples$beta, lambda_mat = samples$lambda, nu_mat = samples$nu,
                                   a_cur = a_vec[k], b_cur = b_vec[k], omega_cur = omega_vec[k], sigmaSq_cur = sigmaSq_vec[k],
                                   a0 = a_vec[resample_iter], b0 = b_vec[resample_iter], omega0 = omega_vec[resample_iter], sigmaSq0 = sigmaSq_vec[resample_iter])
          if (null.a) {
            a_vec = c(a_vec, M.step.IS(weights_vec, samples$nu*a_vec[resample_iter]))
            } else {
              a_vec = c(a_vec, a)
            }
          
          if (null.b) {
            b_vec = c(b_vec, M.step.IS(weights_vec, samples$lambda*b_vec[resample_iter]))
            } else {
             b_vec = c(b_vec, b)
            }
          
          if (null.sigmaSq) {
            sigmaSq_vec = c(sigmaSq_vec, M.step_sigmaSq.IS(weights_vec, samples$beta, X, y, n, diagX = diagX))
            } else {
              sigmaSq_vec = c(sigmaSq_vec, sigmaSq)
            }
          
          if (null.omega) {
            omega = M.step_omega.IS(weights_vec, samples$beta, samples$lambda, samples$nu)
            omega_vec = c(omega_vec, omega)
            } else {
              omega_vec = c(omega_vec, omega)
            }
        }
        
        ### convergence?
        if (converge_rule == "hyper") {
          if (isConverged_hyper(a_vec, b_vec, sigmaSq_vec, omega_vec, delta1, delta2, window_size)) {
            cat("Convergence achieved after", k, "iterations.\n")
            break
            }
          } else {
            if (isConverged_rho(a_vec, b_vec, omega_vec, delta3, window_size)) {
              cat("Convergence achieved after", k, "iterations.\n")
              break
            }
          }
        cat(k, "-th iteration completed.\n")
        cat("a=", a_vec[k+1], "\n")
        cat("b=", b_vec[k+1], "\n")
        cat("omega=", omega_vec[k+1], "\n")
        cat("sigmaSq=", sigmaSq_vec[k+1], "\n")
      }
    }
    
    ### produce posterior samples after convergence
    l = length(a_vec)
    a.est = a_vec[l]
    b.est = b_vec[l]
    sigmaSq.est = sigmaSq_vec[l]
    omega.est = omega_vec[l]
    ### produce samples after "stationary"
    output = getsamples.emp(
      num = output_samples, X = X, y = y,
      a = a.est, b = b.est, omega = omega.est, sigmaSq = sigmaSq.est,
      beta0 = beta0, nu0 = nu0, lambda0 = lambda0, burn = output_burnin,
      xi0 = xi0, woodbury = woodbury, approx = approx, diagX = diagX
      )

    log_lik_mat = log_lik(X = X, y = y, beta_samples = output$beta, sigmaSq = sigmaSq.est)

    return(list(a = a_vec, b = b_vec,
                omega = omega_vec, sigmaSq = sigmaSq_vec,
                beta = output$beta,
                log_lik_mat = log_lik_mat))
  }
}


