#' Draw posterior beta under the nu/lambda TPB augmentation
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param phi Variance multiplier in D = phi * nu / lambda. For unit-rate
#'   local scales this is the TPB phi; for mean-one scales use phi * a / b.
#' @param omega Alias for the variance multiplier in the EB parameterization.
#' @param nu,lambda p*1 local-scale vectors
#' @param sigmaSq error variance (scalar)
#' @param local_var Optional local prior variances, used after marginal MH.
#' @param log_local_var Optional log prior variances, preserving very small scales.
#' @param return_log Return beta together with its sign and log absolute value.
#' @param woodbury binary logic variable, apply Woodbury identity or not
#' @param diagX binary logic variable, if TRUE assumes X is diagonal
#' @return Posterior beta vector, or a list with beta, log_abs_beta and sign_beta
#'   when return_log is TRUE.
#' @export
Gibbs_beta <- function(X, y, phi = NULL, sigmaSq = NULL,
                       nu = NULL, lambda = NULL, woodbury = FALSE,
                       diagX = FALSE, omega = NULL, local_var = NULL,
                       log_local_var = NULL, return_log = FALSE) {
  if (is.null(log_local_var)) {
    if (is.null(local_var)) {
      if (is.null(phi)) phi <- omega
      if (is.null(phi) || is.null(nu) || is.null(lambda)) {
        stop("Gibbs_beta requires local_var or phi (or omega), nu, and lambda.")
      }
      local_var <- phi * pmax(nu, .Machine$double.xmin) /
        pmax(lambda, .Machine$double.xmin)
    }
    local_var <- pmax(local_var, .Machine$double.xmin)
    sqrt_d <- sqrt(local_var)
    log_local_var <- log(local_var)
  } else {
    sqrt_d <- exp(log_local_var/2)
  }
  if (is.null(sigmaSq)) stop("Gibbs_beta requires sigmaSq.")
  sigmaSq <- max(sigmaSq, .Machine$double.xmin)
  finish <- function(z) {
    beta <- sqrt_d * as.vector(z)
    if (!return_log) return(beta)
    list(beta=beta,log_abs_beta=log_local_var/2+log(abs(z)),sign_beta=sign(z))
  }

  if (!diagX) {
    n <- nrow(X)
    p <- ncol(X)
    A <- sweep(X / sqrt(sigmaSq), 2, sqrt_d, FUN = "*")
    if (any(!is.finite(A)))
      return(tpb_gaussian_qr(X, y, sigmaSq, log_local_var, return_log))
    if (!woodbury) {
      gram <- crossprod(A)
      loses_unit <- any(diag(gram) + 1 == diag(gram))
      R <- if (loses_unit) NULL else
        tryCatch(chol(gram + diag(p)), error = function(e) NULL)
      if (is.null(R)) return(tpb_gaussian_qr(X, y, sigmaSq, log_local_var, return_log))
      mean_z <- backsolve(R, forwardsolve(t(R), crossprod(A, y / sqrt(sigmaSq))))
      beta <- finish(as.vector(mean_z + backsolve(R, rnorm(p))))
    } else {
      u <- rnorm(p)
      delta <- rnorm(n)
      rhs <- y / sqrt(sigmaSq) - A %*% u - delta
      gram <- tcrossprod(A)
      loses_unit <- any(diag(gram) + 1 == diag(gram))
      R <- if (loses_unit) NULL else
        tryCatch(chol(gram + diag(n)), error = function(e) NULL)
      if (is.null(R)) {
        # Solve in coefficient coordinates when large scales erase the noise.
        # QR avoids the cancellation in the Woodbury u + correction draw.
        return(tpb_gaussian_qr(X, y, sigmaSq, log_local_var, return_log))
      } else {
        w <- backsolve(R, forwardsolve(t(R), rhs))
      }
      beta <- finish(as.vector(u + crossprod(A, w)))
    }
    return(beta)
  }

  log_s <- pmin(log_local_var/2, .5*log(sigmaSq)-log(abs(diag(X))))
  A <- diag(X)*exp(log_s)/sqrt(sigmaSq)
  precision <- A^2 + exp(2*log_s-log_local_var)
  z <- A*y/sqrt(sigmaSq)/precision + rnorm(length(y))/sqrt(precision)
  beta <- exp(log_s)*z
  if (!return_log) return(beta)
  list(beta=beta, log_abs_beta=log_s+log(abs(z)), sign_beta=sign(z))
}

# QR of the augmented Gaussian precision, without forming X'X.
tpb_gaussian_qr <- function(X, y, sigmaSq, log_var, return_log) {
  p <- ncol(X)
  log_s <- pmin(log_var/2, .5*(log(sigmaSq)-log(colSums(X^2))))
  A <- sweep(X/sqrt(sigmaSq), 2, exp(log_s), "*")
  augmented <- rbind(A, diag(exp(log_s-log_var/2), p))
  decomp <- qr(augmented, LAPACK=TRUE)
  rhs <- qr.qty(decomp, c(y/sqrt(sigmaSq), rep(0, p)))[seq_len(p)]
  z <- numeric(p)
  z[decomp$pivot] <- backsolve(qr.R(decomp), rhs+rnorm(p))
  beta <- exp(log_s)*z
  if (!return_log) return(beta)
  list(beta=beta, log_abs_beta=log_s+log(abs(z)), sign_beta=sign(z))
}
