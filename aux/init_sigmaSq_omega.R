#' Estimate Omega using Method of Moments given estimate of sigmaSq

#' @param X Numeric matrix. The design matrix of dimension n x p.
#' @param y Numeric vector. The response vector of length n.
#' @param sigmaSq_hat Numeric scalar. The estimate of error variance (e.g., from ridge or organic lasso).
#' @param epsilon Numeric scalar. A small positive value to return if the estimate is non-positive. Default is 1e-6.
#'
#' @return Numeric scalar. The estimated value of omega.
#' @export
estimate_omega_mom <- function(X, y, sigmaSq_hat, epsilon = 1e-6) {
  n <- nrow(X)
  y_sq_sum <- sum(y^2)
  trace_XX <- sum(X^2)
  
  # From: 1/n * y'y = 1/n * trace(X'X) * omega + sigmaSq
  # omega = (y'y - n * sigmaSq) / trace(X'X)
  numerator <- y_sq_sum - (n * sigmaSq_hat)
  
  # Check for negative variance estimate
  # This can happen if sigmaSq_hat is overestimated or signal is very weak.
  if (numerator <= 0) {
    warning("Method of Moments estimate for omega is non-positive. Returning epsilon.")
    return(epsilon)
  }

  omega_hat <- numerator / trace_XX
  
  return(omega_hat)
}

#' Initialization Wrapper
#' A wrapper to bundle the sigmaSq and omega initialization.
#' 
#' @param X Numeric matrix.
#' @param y Numeric vector.
#' @param option Binary variable, either "olasso" for oracle lasso with lam = log(p)/n or "ridge" based on ridge regression with CV.
#' 
#' @return A list containing initialized omega and sigmaSq.
#' @importFrom natural olasso
#' @importFrom glmnet cv.glmnet
#' @export
get_initial_values <- function(X, y, option = "ridge") {
  
  if (option == "ridge") {
    cv_ridge <- glmnet::cv.glmnet(X, y, alpha = 0, intercept = FALSE)
    best_lambda <- cv_ridge$lambda.min
    # Get coefficients and residuals to estimate sigmaSq
    beta_ridge <- as.numeric(coef(cv_ridge, s = "lambda.min")[-1])
    y_hat_ridge <- predict(cv_ridge, s = "lambda.min", newx = X)
    sigmaSq_hat <- max(1e-6, var(y - y_hat_ridge))
  } else if (option == "olasso"){
    olasso <- natural::olasso(X, y, intercept = FALSE)
    sigmaSq_hat <- olasso$sig_obj_1
  } else {
    stop("Invalid option for initialization. Please choose 'ridge' or 'olasso'.")
  }
  # Estimate omega using the trace-based Method of Moments
  omega_init <- estimate_omega_mom(X, y, sigmaSq_hat)
  
  return(list(
    sigmaSq = sigmaSq_hat,
    omega = omega_init
  ))
}