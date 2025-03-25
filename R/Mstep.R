#' E-M steps
#' @param sample_matrix a matrix with dimension num_sample*p - nu for estimating a and lambda for estimating b
#' @export
M.step = function(sample_matrix){
  # sample_matrix: a num_sample*p matrix
  empMean_log = colMeans(log(sample_matrix))
  empMean = colMeans(sample_matrix)
  emp_res = mean(empMean - empMean_log) # empirical estimate of 1/p*{\sum E[nu_i] - \sum E[log(nu_i)]}
  f = function(x, emp_res) {
    log(x) + 1 - digamma(x) - emp_res
  }
  #fprime = function(x) {
    #1/x - trigamma(x)
  #}
  solve_for_a = function(emp_res, lower = 1e-3, upper = 100) {
    out = uniroot(f, interval = c(lower, upper), emp_res = emp_res)
    return(out$root)
  }
  solve_for_a(emp_res)
}


#' @param sigmaSq_vec vector, posterior sigmaSq samples
#' @export
M.step_omega = function(beta_matrix, lambda_matrix, nu_matrix, sigmaSq_vec){
  empMean = colMeans(beta_matrix^2 * lambda_matrix/ (nu_matrix * sigmaSq_vec))
  return(mean(empMean))
}
