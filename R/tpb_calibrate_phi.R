#' Calibrate the half-Cauchy scale from a prior effective signal count
#'
#' @param X Design matrix used by the sampler; columns must be nonzero.
#' @param m0 Prior mean effective signal count, strictly between 0 and ncol(X).
#'   For a proportion rho, use m0 = rho * ncol(X).
#' @param hyper_params Gamma hyperprior parameters s_a, r_a, s_b, r_b.
#' @param sigmaSq Optional reference noise variance for calibration. NULL
#'   integrates over the sampler's InvGamma(1.5, 0.5) noise prior.
#' @param n_points Number of deterministic Sobol integration points.
#'   Use a power of two; increase it to check integration accuracy.
#' @details Calibrates phi itself, not sqrt(phi). The target is the prior mean
#'   of sum(d_j * phi * v_j/l_j / (1 + d_j * phi * v_j/l_j)), where
#'   v_j ~ Gamma(a, 1), l_j ~ Gamma(b, 1), and d_j = sum(X[,j]^2)/sigmaSq.
#'   With correlated predictors this is a
#'   coordinate reference count, not prediction degrees of freedom.
#'
#'   Both a and v are integrated analytically. Sobol integration averages over
#'   b, l, the standard half-Cauchy multiplier, an exponential integration
#'   variable, the noise prior (unless sigmaSq is supplied), and column norms.
#'   No response data or MCMC draws are used, and the R random seed is unchanged.
#' @return A list with scale_phi, m0, expected_m0, n_points, sigmaSq, and
#'   integration_change: the absolute change in expected_m0 when using half
#'   the integration points at the same scale. This is a numerical diagnostic,
#'   not an error bound.
#' @export
tpb_calibrate_phi <- function(X, m0,
                              hyper_params = list(s_a = 2.5, r_a = 5,
                                                  s_b = 2.5, r_b = 2.5),
                              sigmaSq = NULL, n_points = 262144) {
  hp <- modifyList(list(s_a = 2.5, r_a = 5, s_b = 2.5, r_b = 2.5), hyper_params)
  if ((!is.null(hp$prior_type_a) && hp$prior_type_a != "gamma") ||
      (!is.null(hp$prior_type_b) && hp$prior_type_b != "gamma"))
    stop("Calibration requires Gamma hyperpriors for a and b.")
  pars <- unlist(hp[c("s_a", "r_a", "s_b", "r_b")])
  if (length(pars) != 4 || any(!is.finite(pars)) || any(pars <= 0))
    stop("Gamma shapes and rates must be positive and finite.")
  norms <- colSums(X^2)
  p <- length(norms)
  if (!p || any(!is.finite(norms)) || any(norms <= 0))
    stop("X must have finite, nonzero columns.")
  if (length(m0) != 1 || !is.finite(m0) || m0 <= 0 || m0 >= p)
    stop("m0 must be strictly between 0 and ncol(X).")
  if (!is.null(sigmaSq) &&
      (length(sigmaSq) != 1 || !is.finite(sigmaSq) || sigmaSq <= 0))
    stop("sigmaSq must be a positive reference variance, or NULL.")
  if (length(n_points) != 1 || !is.finite(n_points) ||
      n_points < 2 || n_points %% 2 != 0)
    stop("n_points must be a positive even integer.")

  u <- gsl::qrng_get(gsl::qrng_alloc("sobol", 7), n_points)
  b <- qgamma(u[,1], shape = hp$s_b, rate = hp$r_b)
  # Gamma(b) = Gamma(b+1) * U^(1/b), kept on the log scale for small b.
  log_lambda <- log(qgamma(u[,2], shape = b+1)) + log(u[,3])/b
  offset <- log(tan(pi*u[,4]/2)) + log(qexp(u[,5])) - log_lambda
  if (is.null(sigmaSq)) {
    offset <- offset + log(qgamma(u[,6], shape = 1.5, rate = 0.5))
  } else {
    offset <- offset - log(sigmaSq)
  }
  # Averaging over columns avoids an n_points by p integration matrix.
  offset <- offset + sort(log(norms))[floor(p*u[,7])+1L]

  contribution <- function(log_scale) {
    ell <- tpb_softplus(log_scale + offset)
    # Integrate v | a and then a: E(exp(-x*v)) = (1 + log(1+x)/r_a)^(-s_a).
    -expm1(-hp$s_a * log1p(ell/hp$r_a))
  }
  mean_count <- function(log_scale) p * mean(contribution(log_scale))
  root <- uniroot(function(z) mean_count(z)-m0, c(-30, 30),
                  extendInt = "upX", tol = 1e-8)$root
  scale_phi <- exp(root)
  if (!is.finite(scale_phi) || scale_phi <= 0)
    stop("Calibrated scale_phi is outside the representable range.")
  values <- contribution(root)
  expected_m0 <- p * mean(values)
  integration_change <- abs(p * mean(values[seq_len(n_points/2)]) - expected_m0)
  if (integration_change > 0.01 * min(m0, p-m0))
    warning("Calibration changes by more than 1% of the smaller effective count; increase n_points.")

  list(scale_phi = scale_phi, m0 = m0, expected_m0 = expected_m0,
       integration_change = integration_change, n_points = n_points,
       sigmaSq = sigmaSq)
}
