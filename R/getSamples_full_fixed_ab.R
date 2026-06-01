#' Fully Bayesian Gibbs Sampler for the TPB Normal Model with Fixed a and b
#' @import GIGrvg
#' @param num_burnin The number of MCMC iterations to discard as burn-in.
#' @param num_samples The number of MCMC iterations to save after the burn-in period.
#' @param fix_phi Logical. If TRUE, phi is fixed to the provided phi_val without sampling.
#' @param phi_val Numeric. The fixed value for phi (if fix_phi = TRUE) or the initial starting value.
#' @param scale_phi Half-Cauchy scale for the typical phi value. Internally this
#'   uses sqrt(phi) ~ C+(0, sqrt(scale_phi)) for a conjugate update.
#' @return list of posterior samples
#' @return list of posterior samples
#' @export

tpb_fullyBayes_fixab <- function(X, y, a, b, num_burnin, num_samples,
                                 woodbury = TRUE, diagX = FALSE,
                                 fix_phi = FALSE, phi_val = NULL,
                                 scale_phi = 1) {
  n <- nrow(X)
  p <- ncol(X)
  num_total <- num_burnin + num_samples
  if (is.null(scale_phi) || !is.finite(scale_phi) || scale_phi <= 0) {
    scale_phi <- 1
  }
  
  beta_samples <- matrix(NA, nrow = num_samples, ncol = p)
  nu_samples <- matrix(NA, nrow = num_samples, ncol = p)
  lambda_samples <- matrix(NA, nrow = num_samples, ncol = p)
  sigmaSq_samples <- numeric(num_samples)
  phi_samples <- numeric(num_samples)
  
  ### Initial values
  sigmaSq_cur <- var(y)
  phi_cur <- if (is.null(phi_val)) scale_phi else phi_val
  # beta_cur <- rnorm(p, 0, sigma = sqrt(sigmaSq_cur))
  nu_cur <- rgamma(p, 1, 1)
  lambda_cur <- rgamma(p, 1, 1)
  
  idx_store <- 0
  
  for (i in 1:num_total) {
    if (!diagX) {
      if (!woodbury) {
        eta_cur <- phi_cur*nu_cur/lambda_cur
        Xrteta <- (X/sqrt(sigmaSq_cur))%*%diag(sqrt(eta_cur))
        AA <- eigen(t(Xrteta) %*% Xrteta)
        nonnegative_eigenVals <- ifelse(AA$values + 1 > 0, 1/(AA$values + 1), 0)
        inverse <- diag(sqrt(eta_cur))%*%AA$vectors%*%diag(nonnegative_eigenVals)%*%t(AA$vectors)%*%diag(sqrt(eta_cur))
        mean_beta <- inverse %*% t(X) %*% (y/sigmaSq_cur)
        beta_cur <- diag(sqrt(eta_cur))%*%AA$vectors%*%diag(sqrt(nonnegative_eigenVals))%*%rnorm(length(mean_beta)) + mean_beta
      } else {
        d_cur <- phi_cur*nu_cur/lambda_cur
        delta <- rnorm(p, 0, sqrt(d_cur))
        # sample auxiliary u ~ N(0, I_n)
        u <- rnorm(n, 0, 1)
        XD <- X %*% delta
        rhs <- (y - XD) / sqrt(sigmaSq_cur) - u
        X_tilde <- sweep(X, 2, sqrt(d_cur), FUN = '*')
        M <- diag(n) + (X_tilde %*% t(X_tilde)) / sigmaSq_cur
        
        eig <- eigen(M, symmetric = TRUE)
        Q <- eig$vectors
        lambdas <- eig$values
        inv_lambdas <- ifelse(lambdas > .Machine$double.eps^0.5, 1/lambdas, 0)
        v <- Q %*% (inv_lambdas * (t(Q) %*% rhs))
        correction_term <- (d_cur * (t(X) %*% v)) / sqrt(sigmaSq_cur)
        beta_cur <- delta + correction_term
      }
    } else {
      mean_beta <- ((diag(X) * y) / sigmaSq_cur)/(diag(X)^2/sigmaSq_cur + lambda_cur/(phi_cur*nu_cur))
      var_beta <- 1/(diag(X)^2/sigmaSq_cur + lambda_cur/(phi_cur*nu_cur))
      beta_cur <- rnorm(length(mean_beta))*sqrt(var_beta) + mean_beta
    }
    
    for (d in 1:p){
      chi <- max(exp(2 * log(abs(beta_cur[d])) + log(lambda_cur[d]) - log(phi_cur)),
                 10^(-306)) # 307 is smallest integer that didn't result in the code breaking here
      nu_cur[d] <- GIGrvg::rgig(1, lambda = a-0.5, chi = chi, psi = 2)
    }
    rate_lam <- 1 + beta_cur^2 / (2*phi_cur*nu_cur)
    lambda_cur <- rgamma(p, shape = b+0.5, rate  = rate_lam)
    
    # sigmaSq | beta,y ~ Inv-Gamma(n/2, SSE/2)
    resid <- y - X %*% beta_cur
    sigmaSq_cur <- 1 / rgamma(1, shape = n/2, rate  = sum(resid^2)/2)
    if (i %% 2000 == 0) {
      print(sum(resid^2)/2)
    }
    if (!fix_phi) {
      w <- 1 / rgamma(1, shape = 1, rate  = 1/scale_phi + 1/phi_cur)
      phi_cur <- 1 / rgamma(1, shape = (p+1)/2, rate = 1/w+sum(lambda_cur * beta_cur^2/(2*nu_cur)))
    }
    
    # store after burn-in
    if(i > num_burnin) {
      idx_store <- idx_store + 1
      beta_samples[idx_store, ] <- beta_cur
      lambda_samples[idx_store, ] <- lambda_cur
      nu_samples[idx_store, ] <- nu_cur
      sigmaSq_samples[idx_store] <- sigmaSq_cur
      phi_samples[idx_store] <- phi_cur
    }
  }
  return(
    list(beta = beta_samples,
         lambda = lambda_samples,
         nu = nu_samples,
         sigmaSq = sigmaSq_samples,
         phi = phi_samples)
  )
}
