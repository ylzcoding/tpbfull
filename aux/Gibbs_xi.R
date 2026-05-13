#' @param p ncol(X)
#' @param a hyper-parameter, real number
#' @param nu p*1 coefficient vector
#' @export
Gibbs_xi = function(a, p, nu){return(rgamma(p, shape = a, rate = 1/nu))}
