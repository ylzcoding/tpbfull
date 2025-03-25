#' define a convergence criterion
#' @param vec1 vector of a
#' @param vec2 vector of b
#' @param vec3 vector of omega
#' @param delta1 numeric, generally fix delta1 at 0.001 to avoid devision by 0
#' @param delta2 numeric
#' @param window_size numeric, size of the sliding window
#' @export
isConverged_hyper = function(vec1, vec2, vec3, delta1, delta2, window_size){
  if(length(vec1) < window_size) {
    return(FALSE)
  }
  # maximum last three consecutive relative errors
  relative_err1 = max(abs(diff(tail(vec1, window_size)))/ (abs(head(tail(vec1, window_size), -1)) + delta1))
  relative_err2 = max(abs(diff(tail(vec2, window_size)))/ (abs(head(tail(vec2, window_size), -1)) + delta1))
  relative_err3 = max(abs(diff(tail(vec3, window_size)))/ (abs(head(tail(vec3, window_size), -1)) + delta1))

  if(max(relative_err1, relative_err2, relative_err3) < delta2){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
