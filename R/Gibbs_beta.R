#' draw posterior beta
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param nu p*1 coefficient vector
#' @param lambda p*1 coefficient vector
#' @param omega coefficient, real number
#' @param sigmasq coefficient, real number
#' @return posterior vector beta
#' @export
Gibbs_beta = function(X, y, omega, sigmaSq, nu, lambda){
  # nu, lambda: p*1 coefficient vector, omega: real number
  # beta_j ~ N(0, sigmaSq*omega*nu_j/lambda_j)
  # reparameterization: omega = phi*a/b
  eta <- omega*nu/lambda
  Xrteta <- X%*%diag(sqrt(eta))
  AA <- eigen(t(Xrteta) %*% Xrteta)
  nonnegative_eigenVals <- ifelse(AA$values + 1 > 0, 1/(AA$values + 1), 0)
  inverse <- diag(sqrt(eta))%*%AA$vectors%*%diag(nonnegative_eigenVals)%*%t(AA$vectors)%*%diag(sqrt(eta))
  mean_beta <- inverse %*% t(X) %*% y

  return(sqrt(sigmaSq)*diag(sqrt(eta))%*%AA$vectors%*%diag(sqrt(nonnegative_eigenVals))%*%rnorm(length(mean_beta)) + mean_beta)
}
