#' Metropolis-Hastings sampler, assuming an uniform prior unif(lower=0, upper)
#' @param vec a vector that will be used, nu for sampling a and lambda for sampling b
#' @param step numeric, step size
#' @return posterior samples and acceptance rates
#' @export
mh = function(p, x0, vec, step, lower = 0, upper = 1){
  x_star = x0 + rnorm(1, 0, step) # new_proposal
  if (x_star > lower & x_star < upper) {
    alpha = exp(sum(dgamma(vec, shape = x_star, rate = x_star, log = TRUE) - dgamma(vec, shape = x0, rate = x0, log = TRUE)))
  } else {
    alpha <- 0
  }
  accept = min(1, alpha)
  U = runif(1, 0, 1)
  gen = ifelse(U<=accept, x_star, x0)
  return(list(gen = gen, acc = accept))
}
