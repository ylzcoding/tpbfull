# Gamma and inverse-Gaussian draws kept on the log scale.
tpb_log_gamma <- function(n, shape) {
  shape <- rep_len(shape, n)
  small <- shape < 1
  draw <- log(rgamma(n, shape + small))
  draw[small] <- draw[small] + log(runif(sum(small))) / shape[small]
  draw
}

tpb_log_add <- function(x, y) pmax(x, y) + log1p(exp(-abs(x - y)))

tpb_log_ig <- function(log_mean, log_shape) {
  n <- length(log_mean)
  log_w <- log_mean + log(rnorm(n)^2) - log(2) - log_shape
  log_den <- numeric(n)
  large <- log_w >= 0
  v <- exp(-log_w[large])
  log_den[large] <- log_w[large] + log(1 + v + sqrt(1 + 2*v))
  w <- exp(log_w[!large])
  log_den[!large] <- log1p(w + sqrt(w*(2+w)))
  lower <- runif(n) <= plogis(log_den)
  log_mean + ifelse(lower, -log_den, log_den)
}

tpb_log_hyperprior <- function(lp, hp) {
  prior_a <- switch(hp$prior_type_a,
    gamma = (hp$s_a-1)*lp[1] - hp$r_a*exp(lp[1]),
    hcauchy = -tpb_softplus(2*(lp[1]-log(hp$scale_a))))
  prior_b <- switch(hp$prior_type_b,
    gamma = (hp$s_b-1)*lp[2] - hp$r_b*exp(lp[2]),
    hcauchy = -tpb_softplus(2*(lp[2]-log(hp$scale_b))))
  prior_a + prior_b - tpb_softplus(2*(lp[3]-log(hp$scale_phi)))
}

tpb_conditional_locals <- function(log_beta, log_lambda, log_xi, a, b, phi) {
  log_omega <- log(phi) + log(a) - log(b)
  log_q <- tpb_log_add(2*log_beta + log_lambda - log_omega, log(2) + log_xi)
  log_nu <- tpb_log_ig((log_q-log(2*a))/2, log_q)
  log_xi <- tpb_log_gamma(length(log_nu), a) + log_nu
  log_rate <- tpb_log_add(log(b), 2*log_beta-log(2)-log_omega-log_nu)
  log_lambda <- tpb_log_gamma(length(log_nu), b+.5) - log_rate
  list(nu = log_nu, lambda = log_lambda, xi = log_xi)
}

tpb_conditional_target <- function(log_beta, log_nu, log_lambda, log_xi, hp) {
  p <- length(log_nu)
  terms <- 2*log_beta + log_lambda - log_nu
  peak <- max(terms)
  log_quad <- peak + log(sum(exp(terms-peak)))
  sum_nu <- sum(exp(log_nu)); sum_lambda <- sum(exp(log_lambda))
  sum_log_lambda <- sum(log_lambda); sum_log_xi <- sum(log_xi)
  function(lp) {
    param <- exp(lp)
    if (any(!is.finite(param)) || any(param <= 0)) return(-Inf)
    a <- param[1]; b <- param[2]
    log_omega <- lp[3]+lp[1]-lp[2]
    gaussian <- -.5*(p*log_omega + exp(log_quad-log_omega))
    # Gamma(nu;a,a) times Gamma(xi;a,1/nu), omitting theta-constant terms.
    nu_xi <- p*a*lp[1]-2*p*lgamma(a)-a*sum_nu+(a-1)*sum_log_xi
    lambda <- p*b*lp[2]-p*lgamma(b)+(b-1)*sum_log_lambda-b*sum_lambda
    gaussian + nu_xi + lambda + tpb_log_hyperprior(lp, hp) + sum(lp)
  }
}
