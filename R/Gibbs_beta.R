#' Draw posterior beta under the nu/lambda TPB augmentation
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param phi Local beta variance scale. This is omega in the EB sampler,
#'   and phi * a / b inside the fully Bayesian sampler.
#' @param nu,lambda p*1 local-scale vectors
#' @param sigmaSq error variance (scalar)
#' @param woodbury binary logic variable, apply Woodbury identity or not
#' @param diagX binary logic variable, if TRUE assumes X is diagonal
#' @return posterior vector beta
#' @export
Gibbs_beta <- function(X, y, a = NULL, b = NULL, phi = NULL, sigmaSq = NULL,
                       nu = NULL, lambda = NULL, woodbury = FALSE,
                       diagX = FALSE, omega = NULL) {
  if (is.null(phi)) {
    phi <- omega
  }
  if (is.null(phi) || is.null(sigmaSq) || is.null(nu) || is.null(lambda)) {
    stop("Gibbs_beta requires phi (or omega), sigmaSq, nu, and lambda.")
  }

  phi <- max(phi, .Machine$double.xmin)
  sigmaSq <- max(sigmaSq, .Machine$double.xmin)
  nu <- pmax(nu, .Machine$double.xmin)
  lambda <- pmax(lambda, .Machine$double.xmin)
  local_var <- pmax(phi * nu / lambda, .Machine$double.xmin)

  if (!diagX) {
    if (!woodbury) {
      eta <- local_var
      Xrteta <- (X / sqrt(sigmaSq)) %*% diag(sqrt(eta))
      AA <- eigen(t(Xrteta) %*% Xrteta)
      nonnegative_eigenVals <- ifelse(AA$values + 1 > 0, 1 / (AA$values + 1), 0)
      inverse <- diag(sqrt(eta)) %*% AA$vectors %*%
        diag(nonnegative_eigenVals) %*% t(AA$vectors) %*% diag(sqrt(eta))
      mean_beta <- inverse %*% t(X) %*% (y / sigmaSq)
      beta <- diag(sqrt(eta)) %*% AA$vectors %*%
        diag(sqrt(nonnegative_eigenVals)) %*% rnorm(length(mean_beta)) + mean_beta
    } else {
      n <- nrow(X)
      p <- ncol(X)
      d <- local_var
      u <- rnorm(p, 0, sqrt(d))
      delta <- rnorm(n, 0, 1)
      XD <- sweep(X, 2, d, FUN = "*")
      M <- XD %*% t(X) + sigmaSq * diag(n)

      eig <- eigen(M, symmetric = TRUE)
      M_eig_vals_inv <- ifelse(eig$values > 1e-8, 1 / eig$values, 0)
      M_inv <- eig$vectors %*% diag(M_eig_vals_inv) %*% t(eig$vectors)
      rhs <- y - X %*% u - sqrt(sigmaSq) * delta
      w <- M_inv %*% rhs
      beta <- u + sweep(t(XD) %*% w, 1, 1, FUN = "*")
    }
    return(beta)
  }

  prior_precision <- 1 / local_var
  mean_beta <- ((diag(X) * y) / sigmaSq) /
    (diag(X)^2 / sigmaSq + prior_precision)
  var_beta <- 1 / (diag(X)^2 / sigmaSq + prior_precision)
  rnorm(length(mean_beta)) * sqrt(var_beta) + mean_beta
}
