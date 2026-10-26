#' Fully Bayesian TPB Sampler with Local eta_j ~ Uniform(0, 1)
#'
#' Baseline sampler inspired by the Yang-style TPB parameterization
#' a_j = eta_j and b_j = 1 - eta_j. This keeps the FPB extra shape parameter
#' fixed at sigma_j = 1, so the sampler remains a three-parameter TPB model
#' with local shape learning rather than the full four-parameter beta prior.
#' Uses sigma^2 ~ InvGamma(1.5, 0.5) and phi ~ half-Cauchy(0, scale_phi),
#' matching the main TPB sampler.
#'
#' @param X n by p design matrix.
#' @param y Response vector of length n.
#' @param num_burnin Number of burn-in iterations.
#' @param num_samples Number of posterior samples to save.
#' @param woodbury Logical; use Woodbury identity in beta updates.
#' @param diagX Logical; treat X as diagonal.
#' @param fix_phi Logical; if TRUE, phi is fixed at phi_val.
#' @param phi_val Fixed or initial value for phi.
#' @param scale_phi Half-Cauchy scale for phi (not sqrt(phi)). Defaults to 1.
#' @param eta_init Initial eta value. Either scalar or vector of length p.
#' @param mh_step_eta Random-walk proposal sd on logit(eta).
#' @param seed Optional random seed.
#' @param sampler "conditional" retains nu, lambda and xi; "marginal" integrates
#'   the local scales out of the eta/phi update and refreshes them exactly.
#' @param interweave Add the transformed, beta-collapsed update for the marginal sampler.
#' @param interweave_every Apply the additional block once every this many iterations.
#' @param adapt_slice Tune phi slice width during burn-in, then freeze it.
#' @param store_locals Save nu and lambda draws for the conditional sampler.
#' @details The additional block uses a joint MH proposal on the local logit(eta)
#'   vector and a slice update of log(phi), then immediately redraws beta.
#' @return Posterior draws for beta, eta, a, b, sigmaSq and phi; acceptance rates
#'   and final slice width. nu and lambda are NULL unless store_locals = TRUE.
#' @export
tpb_fullyBayes_eta_uniform <- function(X, y,
                                       num_burnin,
                                       num_samples,
                                       woodbury = TRUE,
                                       diagX = FALSE,
                                       fix_phi = FALSE,
                                       phi_val = NULL,
                                       scale_phi = 1,
                                       eta_init = 0.5,
                                       mh_step_eta = 2.5,
                                       seed = NULL,
                                       sampler = c("conditional", "marginal"),
                                       interweave = sampler[1] == "marginal",
                                       interweave_every = 1,
                                       adapt_slice = FALSE,
                                       store_locals = sampler[1] == "conditional") {
  if (!is.null(seed)) set.seed(seed)
  sampler <- match.arg(sampler)
  if (sampler == "conditional" && interweave)
    stop("Use interweave = FALSE with the conditional sampler.")
  if (sampler == "marginal" && store_locals)
    stop("Use store_locals = FALSE with the marginal sampler.")
  if (length(interweave_every) != 1 || !is.finite(interweave_every) ||
      interweave_every < 1 || interweave_every != floor(interweave_every))
    stop("interweave_every must be a positive integer.")
  n <- nrow(X); p <- ncol(X)
  num_total <- num_burnin + num_samples
  eta_cur <- rep_len(eta_init, p)
  if (!(length(eta_init) %in% c(1, p)) || any(!is.finite(eta_cur)) ||
      any(eta_cur <= 0 | eta_cur >= 1)) stop("eta_init must lie in (0, 1).")
  if (!is.finite(mh_step_eta) || mh_step_eta <= 0) stop("mh_step_eta must be positive.")
  if (!is.finite(scale_phi) || scale_phi <= 0) stop("scale_phi must be positive.")
  phi_cur <- if (is.null(phi_val)) scale_phi else phi_val
  if (!is.finite(phi_cur) || phi_cur <= 0) stop("phi_val must be positive.")
  h <- qlogis(eta_cur)
  shape <- tpb_eta_shapes(h)
  sigmaSq_cur <- 1 / rgamma(1, shape = 1.5, rate = .5)
  log_nu <- log(rgamma(p, 1)); log_lambda <- log(rgamma(p, 1))
  log_var <- log(phi_cur) + log_nu - log_lambda
  if (sampler == "conditional") log_xi <- tpb_log_gamma(p, shape$a) + log_nu

  beta_samples <- eta_samples <- b_samples <- matrix(NA_real_, num_samples, p)
  nu_samples <- lambda_samples <- if (store_locals) matrix(NA_real_, num_samples, p) else NULL
  sigmaSq_samples <- phi_samples <- numeric(num_samples)
  eta_accept <- numeric(p)
  block_accept <- block_count <- 0L
  slice_width <- 1
  slice_jumps <- 0
  slice_updates <- 0L

  for (iter in seq_len(num_total)) {
    gaussian <- Gibbs_beta(X, y, sigmaSq = sigmaSq_cur, log_local_var = log_var,
                           return_log = TRUE, woodbury = woodbury, diagX = diagX)
    beta_cur <- gaussian$beta; log_beta <- gaussian$log_abs_beta
    sigmaSq_cur <- Gibbs_sigmaSq(n, X, y, beta_cur)
    proposal <- h + rnorm(p, sd = mh_step_eta)
    proposed_shape <- tpb_eta_shapes(proposal)
    log_phi_before <- NULL

    if (sampler == "conditional") {
      log_q <- tpb_log_add(2*log_beta + log_lambda - log(phi_cur), log(2) + log_xi)
      log_nu <- tpb_log_ig((log_q-log(2))/2, log_q)
      log_xi <- tpb_log_gamma(p, shape$a) + log_nu
      log_lambda <- tpb_log_gamma(p, shape$b+.5) -
        tpb_softplus(2*log_beta-log(2)-log(phi_cur)-log_nu)
      target_eta <- function(shape) {
        # b - 1 = -a; the latter avoids cancellation when eta is tiny.
        (shape$a-1)*log_xi - 2*lgamma(shape$a) - shape$a*log_lambda -
          lgamma(shape$b) + log(shape$a) + log(shape$b)
      }
      log_acc <- target_eta(proposed_shape) - target_eta(shape)
    } else {
      target_eta <- function(shape) {
        tpb_log_marginal_beta(beta_cur, shape$a, shape$b, phi_cur, log_beta) +
          log(shape$a) + log(shape$b)
      }
      log_acc <- target_eta(proposed_shape) - target_eta(shape)
    }
    accept <- !is.na(log_acc) & log(runif(p)) < log_acc
    h[accept] <- proposal[accept]
    eta_accept <- eta_accept + accept
    shape <- tpb_eta_shapes(h)

    if (sampler == "conditional") {
      if (!fix_phi) {
        log_phi_before <- log(phi_cur)
        target_phi <- function(lp) {
          if (!is.finite(exp(lp)) || exp(lp) <= 0) return(-Inf)
          -.5*(p*lp + sum(exp(2*log_beta+log_lambda-log_nu-lp))) -
            tpb_softplus(2*(lp-log(scale_phi))) + lp
        }
        phi_cur <- exp(tpb_slice(log_phi_before, target_phi, slice_width))
      }
      log_var <- log(phi_cur) + log_nu - log_lambda
    } else {
      if (!fix_phi) {
        target_phi <- function(lp) {
          sum(tpb_log_marginal_beta(beta_cur, shape$a, shape$b, exp(lp), log_beta)) -
            tpb_softplus(2*(lp-log(scale_phi))) + lp
        }
        log_phi <- log(phi_cur)
        proposal_phi <- rnorm(1, log_phi, .1)
        log_acc <- target_phi(proposal_phi) - target_phi(log_phi)
        if (!is.na(log_acc) && log(runif(1)) < log_acc) phi_cur <- exp(proposal_phi)
      }
      log_var <- tpb_draw_local_var(beta_cur, shape$a, shape$b, phi_cur, log_beta, TRUE)
      if (interweave && iter %% interweave_every == 0) {
        log_phi_before <- if (!fix_phi) log(phi_cur) else NULL
        moved <- tpb_eta_interweave(X, y, log_var, h, phi_cur, sigmaSq_cur, scale_phi,
                                    mh_step_eta, slice_width, fix_phi, woodbury, diagX)
        h <- moved$logit_eta; phi_cur <- moved$phi; log_var <- moved$log_var
        beta_cur <- moved$beta
        shape <- tpb_eta_shapes(h)
        block_accept <- block_accept + moved$accepted
        block_count <- block_count + 1L
      }
    }
    if (adapt_slice && iter <= num_burnin && !is.null(log_phi_before)) {
      slice_updates <- slice_updates + 1L
      slice_jumps <- slice_jumps + abs(log(phi_cur)-log_phi_before)
      if (slice_updates %% 100L == 0) {
        if (slice_jumps > 0) slice_width <- 2*slice_jumps/100
        slice_jumps <- 0
      }
    }
    if (iter > num_burnin) {
      k <- iter-num_burnin
      beta_samples[k, ] <- beta_cur
      eta_samples[k, ] <- shape$a; b_samples[k, ] <- shape$b
      sigmaSq_samples[k] <- sigmaSq_cur; phi_samples[k] <- phi_cur
      if (store_locals) {
        nu_samples[k, ] <- exp(log_nu); lambda_samples[k, ] <- exp(log_lambda)
      }
    }
  }
  list(beta = beta_samples, eta = eta_samples, a = eta_samples, b = b_samples,
       lambda = lambda_samples, nu = nu_samples, sigmaSq = sigmaSq_samples,
       phi = phi_samples, eta_acceptance = eta_accept/num_total,
       interweave_acceptance = if (block_count > 0) block_accept/block_count else NA_real_,
       sigma_i = rep(1, p), mh_step_eta = mh_step_eta, sampler = sampler,
       interweave_every = if (interweave) interweave_every else NULL,
       slice_width = slice_width)
}
