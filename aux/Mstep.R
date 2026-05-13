#' E-M steps
#' @import bzinb
#' @param sample_matrix a matrix with dimension num_sample*p - nu for estimating a and lambda for estimating b
#' @param param_current The current value of the parameter being updated (i.e., a_k or b_k)
#' @return The updated value for the parameter (a_{k+1} or b_{k+1}).
#' @export
M.step = function(sample_matrix, param_current){
  # Add a small epsilon for numerical stability before taking the log
  empMean_log = colMeans(log(sample_matrix + .Machine$double.eps))
  # Calculate the mean of these log-means
  mean_of_empMean_log <- mean(empMean_log)
  idigamma_input <- log(param_current) + mean_of_empMean_log
  # The inverse digamma function can fail for very large negative inputs.
  # Add a safeguard.
  if (!is.finite(idigamma_input)) {
    warning("M.step for a/b received non-finite input.")
    # Return a reasonable small value instead of failing.
    return(1e-4) 
  }
  param_new <- bzinb::idigamma(idigamma_input)
  return(max(1e-6, param_new))
}


#' M-step for sigmaSq using a numerically stable and memory-efficient approach.
#'
#' This version directly computes the expectation of the Sum of Squared Residuals,
#' which avoids catastrophic cancellation and guarantees a non-negative result without
#' creating large p*p matrices.
#'
#' @param beta_matrix S*p matrix of posterior beta samples from the E-step.
#' @param X, y, n The data and its dimensions.
#' @param diagX Logical, if TRUE, assumes X is a diagonal matrix.
#' @return A single, non-negative numeric value for the updated sigmaSq.
#' @export
M.step_sigmaSq <- function(beta_matrix, X, y, n, diagX = FALSE) {
  
  num_samples <- nrow(beta_matrix)
  if (!diagX) {
    # For a general X matrix
    # Calculate the sum of squared residuals for each posterior sample of beta.
    # sapply will return a vector of SSR values, one for each sample.
    ssr_per_sample <- sapply(1:num_samples, function(i) {
      # Get the i-th row (i.e., i-th posterior sample for beta)
      beta_i <- beta_matrix[i, ]
      # Calculate the residuals for this beta sample
      residuals_i <- y - X %*% beta_i
      # Return the sum of squares of these residuals
      sum(residuals_i^2)
    })
    
  } else {
    # For a diagonal X matrix, the calculation is much simpler
    # y_i is predicted by d_i * beta_i, where d_i is the i-th diagonal element of X.
    d <- diag(X) # Extract the diagonal elements
    ssr_per_sample <- sapply(1:num_samples, function(i) {
      beta_i <- beta_matrix[i, ]
      residuals_i <- y - d * beta_i
      sum(residuals_i^2)
    })
  }
  
  # The expectation of the SSR is simply the mean over the posterior samples.
  expected_ssr <- mean(ssr_per_sample)
  # The updated sigmaSq is the average expected squared error.
  # This is guaranteed to be non-negative.
  new_sigmaSq <- expected_ssr / n
  return(max(new_sigmaSq, 1e-16))
}




#' @export
M.step_phi = function(beta_matrix, lambda_matrix, nu_matrix){
  # Add a small epsilon to the denominator for numerical stability
  empMean = colMeans(beta_matrix^2 * lambda_matrix / (nu_matrix + .Machine$double.eps))
  # Safeguard against non-finite results before taking the final mean
  if (any(!is.finite(empMean))) {
    warning("Non-finite values detected in M.step_phi. Taking mean of finite values only.")
    empMean <- empMean[is.finite(empMean)]
    if (length(empMean) == 0) return(1e-4) # Return a default if all are non-finite
  }
  return(mean(empMean))
}
