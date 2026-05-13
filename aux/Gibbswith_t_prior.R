#' @export
StuT.gibbs = function(X, y, warmup, num_samples, df = 1){
  n = dim(X)[1]
  p = dim(X)[2]
  beta_samples = matrix(0, num_samples, p)
  sigmaSq_samples = numeric(num_samples)
  ga_samples = matrix(0, num_samples, p)

  sigmaSq = var(y)
  beta = t(rmvnorm(1, sigma = sigmaSq * diag(p)))
  ga = rep(NA, p)
  for (i in 1:(warmup + num_samples)){
    for (k in 1:p){
      ga[k] = 1/rgamma(1, shape = (df+1)/2, rate = (beta[k]^2 + df)/2)
    }
    var_beta = solve(t(X) %*% X / sigmaSq + diag(1/ga))
    mean_beta = var_beta %*% (t(X) %*% y / sigmaSq)
    beta = t(rmvnorm(1, mean_beta, var_beta))
    sigmaSq = 1/rgamma(1, shape = n/2, rate = t(y - X%*%beta)%*%(y - X%*%beta)/2)

    if (i > warmup){
      beta_samples[i - warmup, ] = beta
      sigmaSq_samples[i - warmup] = sigmaSq
      ga_samples[i - warmup, ] = ga
    }
    if (i %% 1000 == 0){cat(i, "samples generated", "\n")}
  }
  return(list(beta = beta_samples, sigmaSq = sigmaSq_samples, ga = ga_samples))
}




