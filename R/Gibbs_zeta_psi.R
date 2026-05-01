#' Sample local auxiliary variable zeta (IGIG Step 1)
#' @param p ncol(X)
#' @param a hyperparameter a
#' @param b hyperparameter b
#' @param psi p*1 vector of local shrinkage ratios
#' @export
Gibbs_zeta <- function(p, a, b, psi) {
  # zeta | psi ~ Inv-Gamma(a + b, 1 + 1/psi)
  rate_zeta <- 1 + 1 / pmax(psi, 1e-16)
  zeta <- 1 / rgamma(p, shape = a + b, rate = rate_zeta)
  bad_idx <- which(!is.finite(zeta) | is.na(zeta))
  if (length(bad_idx) > 0) zeta[bad_idx] <- 1e-16
  
  return(zeta)
}

#' Sample local shrinkage ratio psi (IGIG Step 2)
#' @param p ncol(X)
#' @param b hyperparameter b
#' @param phi global shrinkage parameter
#' @param beta p*1 coefficient vector
#' @param zeta p*1 vector of local auxiliary variables
#' @export
Gibbs_psi <- function(p, b, phi, beta, zeta) {
  # m_i = beta_i^2 / (2 * phi)
  m <- beta^2 / (2 * pmax(phi, 1e-16))
  # psi | zeta, beta ~ Inv-Gamma(b + 0.5, m + 1/zeta)
  rate_psi <- m + 1 / pmax(zeta, 1e-16)
  psi <- 1 / rgamma(p, shape = b + 0.5, rate = rate_psi)
  bad_idx <- which(!is.finite(psi) | is.na(psi))
  if (length(bad_idx) > 0) psi[bad_idx] <- 1e-16
  
  return(psi)
}