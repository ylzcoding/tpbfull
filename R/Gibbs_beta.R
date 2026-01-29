#' draw posterior beta (a, b, phi parameterization)
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param phi global shrinkage parameter (scalar)
#' @param nu p*1 coefficient vector
#' @param lambda p*1 coefficient vector
#' @param sigmaSq error variance (scalar)
#' @param woodbury binary logic variable, apply woodbury identity or not
#' @param diagX binary logic variable, if TRUE assumes X is diagonal
#' @return posterior vector beta
#' @export
Gibbs_beta = function(X, y, a, b, phi, sigmaSq, nu, lambda, woodbury = FALSE, diagX = FALSE){
  # phi: scalar (global variance component)
  # psi: p*1 vector (local variance component)
  # Prior: beta_j ~ N(0, phi * psi_j)

  if (!diagX) {
    if (!woodbury) {
      # Standard Logic (p < n)
      # eta corresponds to the prior variance of beta
      eta <- phi * nu / lambda
      Xrteta <- (X / sqrt(sigmaSq)) %*% diag(sqrt(eta))
      AA <- eigen(t(Xrteta) %*% Xrteta)
      nonnegative_eigenVals <- ifelse(AA$values + 1 > 0, 1/(AA$values + 1), 0)
      inverse <- diag(sqrt(eta)) %*% AA$vectors %*% diag(nonnegative_eigenVals) %*% t(AA$vectors) %*% diag(sqrt(eta))
      mean_beta <- inverse %*% t(X) %*% (y / sigmaSq)
      beta <- diag(sqrt(eta)) %*% AA$vectors %*% diag(sqrt(nonnegative_eigenVals)) %*% rnorm(length(mean_beta)) + mean_beta
    } else {
      # Woodbury Identity Logic (p >> n)
      n <- nrow(X)
      p <- ncol(X)
      d <- phi * nu / lambda
      # sample u ~ N(0, D) where D = diag(d)
      u <- rnorm(p, 0, sqrt(d))
      # sample delta ~ N(0, I_n)
      delta <- rnorm(n, 0, 1)
      # M = XDX' + sigmaSq*I_n
      XD <- sweep(X, 2, d, FUN="*")   # efficiently compute XD = X %*% D by column-multiplication
      M  <- XD %*% t(X) + sigmaSq * diag(n)

      ############ Eigenvalue decomposition ##########
      # solve w = M^{-1}(y - Xu - sigma*delta)
      eig <- eigen(M, symmetric = TRUE)
      M_eig_vals_inv <- ifelse(eig$values > 1e-8, 1/eig$values, 0)
      M_inv <- eig$vectors %*% diag(M_eig_vals_inv) %*% t(eig$vectors)
      rhs <- y - X %*% u - sqrt(sigmaSq) * delta
      w <- M_inv %*% rhs
      beta <- u + sweep(t(XD) %*% w, 1, 1, FUN="*")
    }
    return(beta)
  } else {
    # Diagonal X Case
    omega <- a*phi/b
    mean_beta <- ((diag(X) * y) / sigmaSq)/(diag(X)^2/sigmaSq + lambda/(omega*nu))
    var_beta <- 1/(diag(X)^2/sigmaSq + lambda/(omega*nu))
    return(rnorm(length(mean_beta))*sqrt(var_beta) + mean_beta)
  }

}
