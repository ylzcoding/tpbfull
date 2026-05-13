#' draw posterior beta
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param nu p*1 coefficient vector
#' @param lambda p*1 coefficient vector
#' @param sigmasq coefficient, real number
#' @param omega hyperparameter, real number
#' @param woodbury binary logic variable, apply woodbury identity or not
#' @param approx binary logic variable, apply thresholding approximation or not 
#' @return posterior vector beta
#' @export
Gibbs_beta = function(X, y, omega, sigmaSq, nu, lambda, 
                      woodbury = FALSE, approx = FALSE, diagX = FALSE){
  # nu, lambda: p*1 coefficient vector, omega: real number
  # beta_j ~ N(0, omega*nu_j/lambda_j)
  if (!diagX) {
    if (!woodbury) {
      eta <- omega*nu/lambda
      Xrteta <- (X/sqrt(sigmaSq))%*%diag(sqrt(eta))
      AA <- eigen(t(Xrteta) %*% Xrteta)
      nonnegative_eigenVals <- ifelse(AA$values + 1 > 0, 1/(AA$values + 1), 0)
      inverse <- diag(sqrt(eta))%*%AA$vectors%*%diag(nonnegative_eigenVals)%*%t(AA$vectors)%*%diag(sqrt(eta))
      mean_beta <- inverse %*% t(X) %*% (y / sigmaSq)
      beta <- diag(sqrt(eta))%*%AA$vectors%*%diag(sqrt(nonnegative_eigenVals))%*%rnorm(length(mean_beta)) + mean_beta
    } else {
      n <- nrow(X)
      p <- ncol(X)
      d <- omega*nu/lambda
      # sample u ~ N(0, D)
      u <- rnorm(p, 0, sqrt(d))
      # sample delta ~ N(0, I_n)
      delta <- rnorm(n, 0, 1)
      # M = XDX' + sigmaSq*I_n
      XD <- sweep(X, 2, d, FUN="*")   # efficiently compute XD = X %*% D by column-multiplication
      M  <- XD %*% t(X) + sigmaSq * diag(n)
      # solve w = M^{-1}(y - Xu - sigma*delta)
      
      ############ Eigenvalue decomposition 
      ############ ISSUE: not stable with the involvement of importance sampling
      eig <- eigen(M, symmetric = TRUE)
      M_eig_vals_inv <- ifelse(eig$values > 1e-8, 1/eig$values, 0)
      M_inv <- eig$vectors %*% diag(M_eig_vals_inv) %*% t(eig$vectors) 
      rhs <- y - X %*% u - sqrt(sigmaSq) * delta
      w <- M_inv %*% rhs
      beta <- u + sweep(t(XD) %*% w, 1, 1, FUN="*")
      
      ############ Cholesky decomposition
      ############ ISSUE: M_reg is occasionally not strict positive definite
      #scale <- max(abs(diag(M)))
      #eps <- max(.Machine$double.eps, scale * .Machine$double.eps * 1e-8)
      #M_reg <- M + diag(eps, n)
      #Lm <- chol(M_reg)
      #rhs <- y - X %*% u - sqrt(sigmaSq) * delta
      #w <- backsolve(Lm, forwardsolve(t(Lm), rhs)) 
      #beta <- u + as.numeric(t(XD) %*% w)
    }
    return(beta) 
  } else {
    mean_beta <- ((diag(X) * y) / sigmaSq)/(diag(X)^2/sigmaSq + lambda/(omega*nu))
    var_beta <- 1/(diag(X)^2/sigmaSq + lambda/(omega*nu))
    return(rnorm(length(mean_beta))*sqrt(var_beta) + mean_beta)
  }
}
