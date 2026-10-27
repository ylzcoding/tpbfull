#' Calibrate the half-Cauchy scale from a prior effective signal count
#'
#' @param X Design matrix used by the sampler; columns must be nonzero.
#' @param m0 Prior mean effective signal count, strictly between 0 and ncol(X).
#'   For a proportion rho, use m0 = rho * ncol(X).
#' @param hyper_params Hyperpriors as in fullGibbs(). For model = "fixed",
#'   supply a and b instead. prior_type_phi selects the global half-Cauchy prior.
#' @param model "tpb" for random a,b, "fixed" for fixed a,b, or "utpb" for
#'   eta ~ Uniform(0,1), a = eta, b = 1-eta.
#' @param sigmaSq Optional reference noise variance for calibration. NULL
#'   integrates over the sampler's InvGamma(1.5, 0.5) noise prior.
#' @param n_points Number of deterministic Sobol integration points.
#'   Use a power of two; increase it to check integration accuracy.
#' @details The default is sqrt(phi) ~ half-Cauchy(0, scale_phi). Set
#'   prior_type_phi = "hcauchy" in hyper_params for the earlier prior on phi.
#'   The target is the prior mean
#'   of sum(d_j * phi * v_j/l_j / (1 + d_j * phi * v_j/l_j)), where
#'   v_j ~ Gamma(a, 1), l_j ~ Gamma(b, 1), and d_j = sum(X[,j]^2)/sigmaSq.
#'   With correlated predictors this is a
#'   coordinate reference count, not prediction degrees of freedom.
#'
#'   v is integrated analytically; a is also integrated analytically when its
#'   hyperprior is Gamma. Sobol integration averages over the remaining prior
#'   variables, an exponential integration variable, and column norms.
#'   No response data or MCMC draws are used, and the R random seed is unchanged.
#' @return A list with scale_phi, prior_type_phi, model, m0, expected_m0,
#'   n_points, sigmaSq, hyper_params, and
#'   integration_change: the absolute change in expected_m0 when using half
#'   the integration points at the same scale. This is a numerical diagnostic,
#'   not an error bound.
#' @export
tpb_calibrate_phi <- function(X, m0,
                              hyper_params = list(s_a = 2.5, r_a = 5,
                                                  s_b = 2.5, r_b = 2.5),
                              sigmaSq = NULL, n_points = 262144,
                              model = c("tpb", "fixed", "utpb")) {
  model <- match.arg(model)
  hp <- modifyList(list(prior_type_a = "gamma", prior_type_b = "gamma",
                        s_a = 2.5, r_a = 5, s_b = 2.5, r_b = 2.5,
                        scale_a = 1, scale_b = 1,
                        prior_type_phi = "sqrt_hcauchy"), hyper_params)
  hp$prior_type_phi <- match.arg(hp$prior_type_phi, c("sqrt_hcauchy", "hcauchy"))
  if (model == "tpb") {
    hp$prior_type_a <- match.arg(hp$prior_type_a, c("gamma", "hcauchy"))
    hp$prior_type_b <- match.arg(hp$prior_type_b, c("gamma", "hcauchy"))
    keys <- c(if (hp$prior_type_a == "gamma") c("s_a", "r_a") else "scale_a",
              if (hp$prior_type_b == "gamma") c("s_b", "r_b") else "scale_b")
  } else {
    keys <- if (model == "fixed") c("a", "b") else character()
  }
  pars <- unlist(hp[keys])
  if (length(pars) != length(keys) || any(!is.finite(pars)) || any(pars <= 0))
    stop("Prior shapes, rates and scales must be positive and finite.")
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

  u <- gsl::qrng_get(gsl::qrng_alloc("sobol", 8), n_points)
  gamma_a <- model == "tpb" && hp$prior_type_a == "gamma"
  if (model == "fixed") {
    a <- hp$a; b <- hp$b
  } else if (model == "utpb") {
    a <- u[,1]; b <- 1-a
  } else {
    b <- if (hp$prior_type_b == "gamma") qgamma(u[,1], hp$s_b, hp$r_b) else
      hp$scale_b * tan(pi*u[,1]/2)
    if (!gamma_a) a <- hp$scale_a * tan(pi*u[,8]/2)
  }
  # Gamma(b) = Gamma(b+1) * U^(1/b), kept on the log scale for small b.
  log_lambda <- log(qgamma(u[,2], shape = b+1)) + log(u[,3])/b
  power <- if (hp$prior_type_phi == "sqrt_hcauchy") 2 else 1
  offset <- power*log(tan(pi*u[,4]/2)) + log(qexp(u[,5])) - log_lambda
  if (is.null(sigmaSq)) {
    offset <- offset + log(qgamma(u[,6], shape = 1.5, rate = 0.5))
  } else {
    offset <- offset - log(sigmaSq)
  }
  # Averaging over columns avoids an n_points by p integration matrix.
  offset <- offset + sort(log(norms))[floor(p*u[,7])+1L]

  contribution <- function(log_scale) {
    ell <- tpb_softplus(power*log_scale + offset)
    # Integrate v | a; integrate a as well when it has a Gamma hyperprior.
    if (gamma_a) -expm1(-hp$s_a * log1p(ell/hp$r_a)) else -expm1(-a*ell)
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

  hp$scale_phi <- scale_phi
  list(scale_phi = scale_phi, prior_type_phi = hp$prior_type_phi, model = model,
       m0 = m0, expected_m0 = expected_m0, hyper_params = hp,
       integration_change = integration_change, n_points = n_points,
       sigmaSq = sigmaSq)
}
