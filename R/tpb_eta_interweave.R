tpb_eta_shapes <- function(logit_eta) {
  list(a = exp(-tpb_softplus(-logit_eta)), b = exp(-tpb_softplus(logit_eta)))
}

tpb_eta_collapsed_target <- function(X, y, log_var, logit_eta, phi, sigmaSq, scale_phi) {
  shape <- tpb_eta_shapes(logit_eta)
  x <- log_var - log(phi)
  negative <- x < 0
  z <- x * ifelse(negative, shape$a, shape$b)
  transform <- function(h, log_phi) {
    shape <- tpb_eta_shapes(h)
    x <- z / ifelse(negative, shape$a, shape$b)
    list(a = shape$a, b = shape$b, x = x, log_var = log_phi + x)
  }
  target <- function(h, log_phi) {
    value <- transform(h, log_phi)
    if (any(!is.finite(value$log_var)) || any(value$a <= 0) || any(value$b <= 0) ||
        !is.finite(exp(log_phi)) || exp(log_phi) <= 0) return(-Inf)
    log_z <- sum(-abs(z) - tpb_softplus(-abs(value$x)) - lbeta(value$a, value$b) -
                   log(ifelse(negative, value$a, value$b)))
    tpb_log_marginal_y(X, y, value$log_var, sigmaSq) + log_z +
      sum(log(value$a) + log(value$b)) + log_phi -
      tpb_softplus(2*(log_phi-log(scale_phi)))
  }
  list(target = target, transform = transform)
}

tpb_eta_interweave <- function(X, y, log_var, logit_eta, phi, sigmaSq, scale_phi,
                               mh_step_eta, slice_width, fix_phi,
                               woodbury, diagX) {
  representation <- tpb_eta_collapsed_target(X, y, log_var, logit_eta, phi,
                                             sigmaSq, scale_phi)
  log_phi <- log(phi)
  # One block proposal avoids p Gaussian factorizations for local eta_j's.
  proposal <- logit_eta + rnorm(length(logit_eta), sd = mh_step_eta/sqrt(length(logit_eta)))
  log_acc <- representation$target(proposal, log_phi) -
    representation$target(logit_eta, log_phi)
  accepted <- !is.na(log_acc) && log(runif(1)) < log_acc
  if (accepted) logit_eta <- proposal
  if (!fix_phi) {
    target <- function(value) representation$target(logit_eta, value)
    log_phi <- tpb_slice(log_phi, target, slice_width)
  }
  value <- representation$transform(logit_eta, log_phi)
  gaussian <- Gibbs_beta(X, y, sigmaSq = sigmaSq, log_local_var = value$log_var,
                         woodbury = woodbury, diagX = diagX, return_log = TRUE)
  list(logit_eta = logit_eta, phi = exp(log_phi), log_var = value$log_var,
       beta = gaussian$beta, log_abs_beta = gaussian$log_abs_beta, accepted = accepted)
}
