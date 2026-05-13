#' Define a convergence criterion. Several numeric vectors, each containing at least window_size + 1 elements.
#' The function calculates the relative change over the last window_size + 1 points for each vector individually, 
#' and checks whether all of them fall below the threshold delta2.
#' @param ..., multiple numerical vectors
#' @param delta1, numeric, generally fix delta1 at 0.001
#' @param delta2, numeric
#' @param window_size numeric, size of the sliding window
#' @export
isConverged_hyper <- function(..., delta1, delta2, window_size = 3) {
  vecs <- list(...)
  if (any(sapply(vecs, length) < window_size + 1)) {
    return(FALSE)
  }
  
  rel_err <- function(vec) {
    w <- tail(vec, window_size + 1)
    diffs <- abs(diff(w))
    denom <- abs(head(w, -1)) + delta1
    rel  <- diffs / denom
    if (all(is.na(rel))) return(Inf)
    max(rel, na.rm = TRUE)
  }
  
  rel_errs <- vapply(vecs, rel_err, numeric(1))
  
  return(max(rel_errs, na.rm = TRUE) < delta2)
}


#' define another convergence criterion
#' @param beta_matrix A num_iteration * p matrix, where each row is the posterior mean estimate of beta at each iteration
#' @param window_size numeric, size of the sliding window
#' @export
isConverged_beta = function(beta_matrix, delta1, delta3, window_size = 3){
  num_iter = nrow(beta_matrix)
  if (num_iter < (window_size+1)) {
    return(FALSE)
  }
  relative_err = numeric(window_size)

  for (i in 1:window_size) {
    # the last (window_size + 1) rows
    mse = mean((beta_matrix[num_iter-window_size+i, ] - beta_matrix[num_iter-window_size+i-1, ])^2)
    mean_squares = mean(beta_matrix[num_iter-window_size+i, ]^2)
    relative_err[i] = mse / (mean_squares + delta1)
  }
  if (all(relative_err < delta3)) {
    return(TRUE)
  }else{
    return(FALSE)
  }
}



#' convergence criterion based on the Hellinger distance of shrinkage factor densities
#' @param a_vec, vector of a
#' @param b_vec, vector of b
#' @param omega_vec, vector of omega
#' @param delta, numeric, error bound
#' @param window_size numeric, size of the sliding window
#' @export
isConverged_rho = function(a_vec, b_vec, omega_vec, delta, window_size = 3){
  if (length(a_vec) < window_size + 1) {
    return(FALSE)
  }

  a_tail <- tail(a_vec, window_size + 1)
  b_tail <- tail(b_vec, window_size + 1)
  omega_tail <- tail(omega_vec, window_size + 1)
  k <- length(a_vec)

  for (i in 1:window_size) {
    # the last (window_size + 1) rows
    d <- hlg_dis(a_cur = a_tail[i+1], b_cur = b_tail[i+1], omega_cur = omega_tail[i+1], 
                 a_prev = a_tail[i], b_prev = b_tail[i], omega_prev = omega_tail[i])
    cat(sprintf("Hellinger distance between iteration %d and iteration %d is %f.\n",
                k - window_size + i - 1, k - window_size + i, d))
    if (is.infinite(d) || d >= delta) {
      return(FALSE)
    }
  }
  return(TRUE)
}





