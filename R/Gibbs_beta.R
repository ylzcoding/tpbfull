#' Draw posterior beta under the nu/lambda TPB augmentation
#' @param X n*p design matrix
#' @param y n*1 response vector
#' @param phi Local beta variance scale. This is omega in the EB sampler,
#'   or phi * a / b in the mean-one nu/lambda representation of marginal TPB.
#' @param nu,lambda p*1 local-scale vectors
#' @param sigmaSq error variance (scalar)
#' @param local_var Optional local prior variances, used after marginal MH.
#' @param woodbury binary logic variable, apply Woodbury identity or not
#' @param diagX binary logic variable, if TRUE assumes X is diagonal
#' @return posterior vector beta
#' @export
Gibbs_beta <- function(X, y, a = NULL, b = NULL, phi = NULL, sigmaSq = NULL,
                       nu = NULL, lambda = NULL, woodbury = FALSE,
                       diagX = FALSE, omega = NULL, local_var = NULL) {
  if (is.null(local_var)) {
    if (is.null(phi)) phi <- omega
    if (is.null(phi) || is.null(nu) || is.null(lambda)) {
      stop("Gibbs_beta requires local_var or phi (or omega), nu, and lambda.")
    }
    local_var <- phi * pmax(nu, .Machine$double.xmin) /
      pmax(lambda, .Machine$double.xmin)
  }
  if (is.null(sigmaSq)) stop("Gibbs_beta requires sigmaSq.")
  sigmaSq <- max(sigmaSq, .Machine$double.xmin)
  local_var <- pmax(local_var, .Machine$double.xmin)

  if (!diagX) {
    n <- nrow(X)
    p <- ncol(X)
    sqrt_d <- sqrt(local_var)
    A <- sweep(X / sqrt(sigmaSq), 2, sqrt_d, FUN = "*")
    if (!woodbury) {
      R <- chol(crossprod(A) + diag(p))
      mean_z <- backsolve(R, forwardsolve(t(R), crossprod(A, y / sqrt(sigmaSq))))
      beta <- sqrt_d * as.vector(mean_z + backsolve(R, rnorm(p)))
    } else {
      u <- rnorm(p)
      delta <- rnorm(n)
      rhs <- y / sqrt(sigmaSq) - A %*% u - delta
      gram <- tcrossprod(A)
      loses_unit <- any(diag(gram) + 1 == diag(gram))
      R <- if (loses_unit) NULL else
        tryCatch(chol(gram + diag(n)), error = function(e) NULL)
      if (is.null(R)) {
        # Large local scales can erase I_n when forming A A' numerically.
        # SVD of A preserves the +1 without clipping eigenvalues or adding jitter.
        decomp <- svd(A)
        correction <- decomp$v %*%
          ((decomp$d / (1 + decomp$d^2)) * crossprod(decomp$u, rhs))
        return(sqrt_d * as.vector(u + correction))
      } else {
        w <- backsolve(R, forwardsolve(t(R), rhs))
      }
      beta <- sqrt_d * as.vector(u + crossprod(A, w))
    }
    return(beta)
  }

  prior_precision <- 1 / local_var
  mean_beta <- ((diag(X) * y) / sigmaSq) /
    (diag(X)^2 / sigmaSq + prior_precision)
  var_beta <- 1 / (diag(X)^2 / sigmaSq + prior_precision)
  rnorm(length(mean_beta)) * sqrt(var_beta) + mean_beta
}
