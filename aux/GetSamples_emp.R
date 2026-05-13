#' @param burn numeric, number of burn-in samples
#' @param num numeric, number of samples to generate after burn-in for each hyperparameter
#' @param xi0 initial value for xi, the auxiliary variable, xi0 = NULL if not needed
#' @return a list of posterior samples generated given a and b
#' @export
getsamples.emp = function(num, X, y, a, b, omega, sigmaSq, beta0, nu0, lambda0,
                          burn, xi0 = NULL, woodbury = FALSE, approx = FALSE, diagX = FALSE){
  aug <- !is.null(xi0)
  p = dim(X)[2]
  n = dim(X)[1]
  # num: number of samples to generate after burn-in
  beta_samples = matrix(NA, nrow = num + burn, ncol = p)
  beta_samples[1, ] = beta0
  nu_samples = matrix(NA, nrow = num + burn, ncol = p)
  nu_samples[1, ] = nu0
  lambda_samples = matrix(NA, nrow = num + burn, ncol = p)
  lambda_samples[1, ] = lambda0
  if (aug) {
    xi_samples = matrix(NA, nrow = num + burn, ncol = p)
    xi_samples[1, ] = xi0
  }
  for (i in 2:(num + burn)){
    beta_samples[i, ] = Gibbs_beta(X = X, y = y, omega = omega, sigmaSq = sigmaSq, nu = nu_samples[i-1, ], lambda = lambda_samples[i-1, ], 
                                   woodbury = woodbury, approx = approx, diagX = diagX)
    if (aug) {
      nu_samples[i, ] = aug_Gibbs_nu(p = p, omega = omega, beta = beta_samples[i, ], lambda = lambda_samples[i-1, ], xi = xi_samples[i-1, ], a = a)
      xi_samples[i, ] = Gibbs_xi(a = a, p = p, nu = nu_samples[i, ])
    } else {
      nu_samples[i, ] = Naive_Gibbs_nu(a = a, p = p, omega = omega, beta = beta_samples[i, ], lambda = lambda_samples[i-1, ])
    }
    lambda_samples[i, ] = Gibbs_lambda(b = b, p = p, omega = omega, beta = beta_samples[i, ], nu = nu_samples[i, ])
  }
  res <- list(beta = beta_samples[(burn + 1):(burn + num), ],
              nu = nu_samples[(burn + 1):(burn + num), ],
              lambda = lambda_samples[(burn + 1):(burn + num), ])
  if (aug) {
    res$xi = xi_samples[(burn + 1):(burn + num), ]
  }
  return(res)
}
