#' @return Hellinger distance between the current shrinkage factor density and the previous one.
#' @export
hlg_dis <- function(a_cur, b_cur, omega_cur, a_prev, b_prev, omega_prev) {
  # calculate normalizing constant
  norm_Z <- function(a, b, omega) {
    phi <- b*omega/a
    f <- function(x) {
      x^(b-1) * (1-x)^(a-1) * (1+(phi-1)*x)^(-a-b)
    }
    Z <- integrate(f, lower = 0, upper = 1, stop.on.error = FALSE)
    return(Z$value)
  }
  
  Z_cur <- norm_Z(a_cur, b_cur, omega_cur)
  Z_prev <- norm_Z(a_prev, b_prev, omega_prev)
  
  phi_cur <- b_cur*omega_cur/a_cur
  phi_prev <- b_prev*omega_prev/a_prev
  g <- function(x) {
    x^((b_cur+b_prev)/2-1) * (1-x)^((a_cur+a_prev)/2-1) * (1+(phi_cur-1)*x)^(-(a_cur+b_cur)/2) * (1+(phi_prev-1)*x)^(-(a_prev+b_prev)/2)
  }
  int <- integrate(g, lower = 0, upper = 1, stop.on.error = FALSE)
  if (is.na(int$value) || !is.finite(int$value)) {
    warning(paste("Main integration in hlg_dis failed, likely due to extreme omega. Message:", int$message), call. = FALSE)
    return(Inf) # Return Inf to signal non-convergence.
  }
  res <- int$value / sqrt(Z_cur * Z_prev)
  # The max(..., 0) handles cases where res might be slightly > 1 due to numerical error
  hlg_dis <- sqrt(max(1-res, 0))
  return(hlg_dis)
}