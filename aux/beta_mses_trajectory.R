#' @param beta_matrix A num_iteration * p matrix, where each row is the posterior mean estimate of beta at each iteration
#' @return A vector of length (num_iteration-1)， where i-th element is the mean squared error between beta_matrix[i, ] and beta_matrix[i+1, ].
#' @export
mse_traj <- function(beta_matrix) {

  num_iterations <- nrow(beta_matrix)
  mse_vec <- numeric(num_iterations - 1)

  for (i in 1:(num_iterations - 1)) {
    mse_vec[i] <- mean((beta_matrix[i+1, ] - beta_matrix[i, ])^2)
  }
  return(mse_vec)
}
