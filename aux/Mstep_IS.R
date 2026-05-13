#' @import bzinb
#' @param sample_matrix a num_sample by p posterior sample matrix, sampled using a0, b0, omega0
#' @param IS_weights a weight vector of length num_sample
#' @export
M.step.IS <- function(IS_weights, sample_matrix){
  empMean.IS <- colMeans(log(sample_matrix) * IS_weights)
  return(bzinb::idigamma(mean(empMean.IS)/mean(IS_weights)))
}


#' @param sigmaSq numeric, estimated sigmaSq from last iteration
#' @param IS_weights a weight vector of length num_sample
#' @export
M.step_omega.IS <- function(IS_weights, beta_matrix, lambda_matrix, nu_matrix){
  
  empMean.IS <- colMeans((beta_matrix^2 * lambda_matrix/ nu_matrix) * IS_weights)
  return(mean(empMean.IS) / mean(IS_weights))
}


#' @param beta_matrix beta matrix with dimension num_sample*p
#' @param IS_weights a weight vector of length num_sample
#' @export
M.step_sigmaSq.IS <- function(IS_weights, beta_matrix, X, y, n, diagX = FALSE){
  beta_samples <- t(beta_matrix)
  num_samples <- dim(beta_samples)[2]
  w <- sum(IS_weights)
  
  if (!diagX) {
    beta_beta_prime <- vector("list", length = num_samples)
    for(i in 1:num_samples){
      beta_beta_prime[[i]] <- IS_weights[i] * (beta_samples[, i] %*% t(beta_samples[, i]))
    }
    empMean_beta_beta_prime <- Reduce("+", beta_beta_prime)/ w
    empMean_beta <- colSums(beta_matrix * IS_weights) / w
    return((sum(diag(t(X) %*% X %*% empMean_beta_beta_prime)) - 2 * t(y) %*% X %*% empMean_beta + t(y) %*% y)/ n)
  } else {
    empMean_beta_beta_prime <- colSums(beta_matrix^2 * IS_weights) / w
    empMean_beta <- colSums(beta_matrix * weights) / w
    return((sum(diag(X)^2*empMean_beta_beta_prime) - 2 * sum(y*diag(X)*empMean_beta) + sum(y^2))/ n) 
  }
}

