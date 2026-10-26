# With v ~ Gamma(a, 1), l ~ Gamma(b, 1), beta | v,l ~ N(0, phi*v/l),
# integrating l gives the following log density of x = log(v), up to a constant.
tpb_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

tpb_local_log_kernel <- function(x, a, b, log_c) {
  (a + b) * x - exp(x) - (b + 0.5) *
    (log_c + tpb_softplus(x - log_c))
}

tpb_local_envelope <- function(a, b, log_c) {
  # Positive root of v^2 + (c + 1/2 - a)*v - (a+b)*c = 0,
  # scaled before evaluation to avoid cancellation and overflow.
  a <- rep_len(a, length(log_c)); b <- rep_len(b, length(log_c))
  k <- a - 0.5
  log_scale <- pmax(log_c, log(abs(k)), (log(a + b) + log_c) / 2)
  B <- exp(log_c - log_scale) - sign(k) * exp(log(abs(k)) - log_scale)
  disc <- sqrt(B^2 + 4 * exp(log(a + b) + log_c - 2 * log_scale))
  mode <- numeric(length(log_c))
  pos <- B >= 0
  mode[pos] <- log(2) + log(a[pos] + b[pos]) + log_c[pos] -
    log_scale[pos] - log(B[pos] + disc[pos])
  mode[!pos] <- log_scale[!pos] + log(disc[!pos] - B[!pos]) - log(2)
  peak <- tpb_local_log_kernel(mode, a, b, log_c)

  left_width <- right_width <- rep(1, length(log_c))
  repeat {
    left <- mode - left_width
    h_left <- tpb_local_log_kernel(left, a, b, log_c)
    expand <- h_left > peak - 1
    if (!any(expand)) break
    left_width[expand] <- 2 * left_width[expand]
  }
  repeat {
    right <- mode + right_width
    h_right <- tpb_local_log_kernel(right, a, b, log_c)
    expand <- h_right > peak - 1
    if (!any(expand)) break
    right_width[expand] <- 2 * right_width[expand]
  }
  slope_left <- a + b - exp(left) - (b + 0.5) * plogis(left - log_c)
  slope_right <- a + b - exp(right) - (b + 0.5) * plogis(right - log_c)
  list(mode = mode, peak = peak, left = left, right = right,
       slope_left = slope_left, slope_right = slope_right,
       cut_left = left + (peak - h_left) / slope_left,
       cut_right = right + (peak - h_right) / slope_right)
}

tpb_draw_local_var <- function(beta, a, b, phi,
                               log_abs_beta = log(abs(beta)), return_log = FALSE) {
  # Exact conditional refresh after marginal MH. One ordinary Gibbs sweep
  # using old local scales would not restore their joint conditional law.
  a <- rep_len(a, length(beta)); b <- rep_len(b, length(beta))
  log_c <- 2 * log_abs_beta - log(2) - log(phi)
  x <- numeric(length(beta))
  zero <- log_abs_beta == -Inf
  if (any(zero)) {
    if (any(a[zero] <= 0.5)) stop("Local conditional at beta = 0 is improper when a <= 1/2.")
    x[zero] <- log(rgamma(sum(zero), shape = a[zero] - 0.5))
  }
  nz <- which(!zero)
  if (length(nz)) {
    env <- tpb_local_envelope(a[nz], b[nz], log_c[nz])
    mass_left <- 1 / env$slope_left
    mass_mid <- env$cut_right - env$cut_left
    mass_right <- -1 / env$slope_right
    pending <- seq_along(nz)
    # The log kernel is strictly concave. Its two tangents and its maximum
    # give an exponential / uniform / exponential rejection envelope.
    while (length(pending)) {
      j <- pending
      piece <- runif(length(j)) * (mass_left[j] + mass_mid[j] + mass_right[j])
      u <- runif(length(j))
      draw <- env$cut_left[j] + u * mass_mid[j]
      log_envelope <- numeric(length(j))
      left <- piece < mass_left[j]
      right <- piece > mass_left[j] + mass_mid[j]
      draw[left] <- env$cut_left[j[left]] + log(u[left]) / env$slope_left[j[left]]
      draw[right] <- env$cut_right[j[right]] + log(u[right]) / env$slope_right[j[right]]
      log_envelope[left | right] <- log(u[left | right])
      accept <- log(runif(length(j))) <=
        tpb_local_log_kernel(draw, a[nz[j]], b[nz[j]], log_c[nz[j]]) - env$peak[j] - log_envelope
      x[nz[j[accept]]] <- draw[accept]
      pending <- j[!accept]
    }
  }
  # l | v,beta ~ Gamma(b+1/2, 1+c/v), so phi*v/l = phi*(v+c)/G.
  # Compute the ratio directly: v and l can both underflow even when it is finite.
  log_sum <- pmax(x, log_c) + log1p(exp(-abs(x - log_c)))
  log_var <- log(phi) + log_sum - log(rgamma(length(beta), shape = b + 0.5))
  if (return_log) log_var else exp(log_var)
}

tpb_log_marginal_beta <- function(beta_vec, a, b, phi,
                                  log_abs_beta = log(abs(beta_vec))) {
  if (!is.finite(phi) || phi <= 0) return(rep(-Inf, length(beta_vec)))
  a <- rep_len(a, length(beta_vec)); b <- rep_len(b, length(beta_vec))
  valid <- is.finite(a) & is.finite(b) & a > 0 & b > 0
  log_c <- 2 * log_abs_beta - log(2) - log(phi)
  log_const <- rep(NA_real_, length(beta_vec))
  log_const[valid] <- lgamma(b[valid] + 0.5) - lgamma(a[valid]) - lgamma(b[valid]) -
    0.5 * (log(2 * pi) + log(phi))
  out <- rep(-Inf, length(beta_vec))
  zero <- log_abs_beta == -Inf
  out[zero & valid] <- Inf
  finite_zero <- zero & valid & a > 0.5
  out[finite_zero] <- log_const[finite_zero] + lgamma(a[finite_zero] - 0.5)

  nz <- which(!zero & valid)
  if (!length(nz)) return(out)
  z <- exp(log_c[nz])
  U <- rep(NaN, length(nz))
  representable <- is.finite(z) & z > 0
  if (any(representable)) {
    U[representable] <- suppressWarnings(gsl::hyperg_U(
      b[nz[representable]] + 0.5, 1.5 - a[nz[representable]], z[representable]
    ))
  }
  ok <- is.finite(U) & U > 0
  out[nz[ok]] <- log_const[nz[ok]] + lgamma(a[nz[ok]] + b[nz[ok]]) + log(U[ok])
  for (j in nz[!ok]) {
    # Evaluate the equivalent Gamma-mixture integral on log(v), scaled at
    # its mode. Neither U nor Gamma(a+b) needs to be representable separately.
    env <- tpb_local_envelope(a[j], b[j], log_c[j])
    width <- (env$right - env$left) / 2
    integrand <- function(t) {
      exp(tpb_local_log_kernel(env$mode + width * t, a[j], b[j], log_c[j]) - env$peak)
    }
    integral <- integrate(integrand, -Inf, 0, rel.tol = 1e-8)$value +
      integrate(integrand, 0, Inf, rel.tol = 1e-8)$value
    out[j] <- log_const[j] + env$peak + log(width) + log(integral)
  }
  out
}
