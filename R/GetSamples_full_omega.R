#' @param num_burnin numeric, number of burn-in samples
#' @param num_output numeric, number of samples to generate after burn-in for each hyperparameter
#' @param xi0 initial value for xi, the auxiliary variable, xi0 = NULL if not needed
#' @param omega_hyper hyper-parameter of the prior for omega
#' @return a list of posterior samples generated given a and b
#' @export
getsamples.full_omega = function(num_output, num_burnin,
                                 X, y, a, b, omega_hyper,
                                 omega0, w0, sigmaSq0, beta0, nu0, lambda0, xi0 = NULL) {
  # prepare posterior samples including posterior omega, given a and b
  # used in E-M within Gibbs when we'd like to empirically estimate a and b, but go with fully Bayesian on omega.
  aug = !is.null(xi0) # if aug is TRUE, go with augmented gibbs sampling for nu
  p = dim(X)[2]
  n = dim(X)[1]
  num_total = num_output + num_burnin

  beta_samples = matrix(NA, nrow = num_total, ncol = p)
  beta_samples[1, ] = beta0
  nu_samples = matrix(NA, nrow = num_total, ncol = p)
  nu_samples[1, ] = nu0
  if (aug) {
    xi_samples = matrix(NA, nrow = num_total, ncol = p)
    xi_samples[1, ] = xi0
  }
  lambda_samples = matrix(NA, nrow = num_total, ncol = p)
  lambda_samples[1, ] = lambda0
  omega_samples = rep(NA, num_total)
  omega_samples[1] = omega0
  w_samples = rep(NA, num_total)
  w_samples[1] = w0
  sigmaSq_samples = rep(NA, num_total)
  sigmaSq_samples[1] = sigmaSq0

  for (i in 2:num_total){
    beta_samples[i, ] = Gibbs_beta(X, y, omega_samples[i-1], sigmaSq_samples[i-1], nu_samples[i-1, ], lambda_samples[i-1, ])
    if (aug) {
      nu_samples[i, ] = aug_Gibbs_nu(a, p, omega_samples[i-1], sigmaSq_samples[i-1], beta_samples[i, ], lambda_samples[i-1, ], xi_samples[i-1, ])
      xi_samples[i, ] = Gibbs_xi(a, p, nu_samples[i, ])
    } else {
      nu_samples[i, ] = Naive_Gibbs_nu(a, p, omega_samples[i-1], sigmaSq_samples[i-1], beta_samples[i, ], lambda_samples[i-1, ])
    }
    lambda_samples[i, ] = Gibbs_lambda(b, p, omega_samples[i-1], sigmaSq_samples[i-1], beta_samples[i, ], nu_samples[i, ])
    omega_samples[i] = Gibbs_omega(p, nu_samples[i, ], lambda_samples[i, ], w_samples[i-1], sigmaSq_samples[i-1], beta_samples[i, ])
    w_samples[i] = Gibbs_w(omega_samples[i], omega_hyper)
    sigmaSq_samples[i] = Gibbs_sigmaSq(n, p, y, X, beta_samples[i, ], omega_samples[i], nu_samples[i, ], lambda_samples[i, ])
  }
  res <- list(beta = beta_samples[(num_burnin + 1):num_total, ],
              nu = nu_samples[(num_burnin + 1):num_total, ],
              lambda = lambda_samples[(num_burnin + 1):num_total, ],
              omega = omega_samples[(num_burnin + 1):num_total],
              w = w_samples[(num_burnin + 1):num_total],
              sigmaSq = sigmaSq_samples[(num_burnin + 1):num_total])
  if (aug) {
    res$xi = xi_samples[(num_burnin + 1):num_total, ]
  }
  return(res)
}
