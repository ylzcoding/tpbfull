Fully Bayesian three-parameter beta regression

The default hyperpriors are `a ~ Gamma(2.5, 5)` and `b ~ Gamma(2.5, 2.5)`
(shape, rate). Set `m0` to the prior mean effective signal count to calibrate
`sqrt(phi) ~ half-Cauchy(0, scale_phi)` before sampling. For a proportion, use
`m0 = rho * ncol(X)`.

```r
fit <- tpb_full_pipeline(X, y, m0 = 0.05 * ncol(X),
                         num_chains = 4, num_iter = 40000, num_warmup = 20000,
                         interweave_every = 3, adapt_slice = TRUE)
fit$hyper_params
fit$calibration
```

This runs four chains with 5,000 burn-in and 5,000 retained draws each.
`auto_find = FALSE` is now the default. The calibration uses X and the priors,
including `sigma^2 ~ InvGamma(1.5, 0.5)`, without using y. It averages over a,
b, phi and the local scales, then solves for one fixed half-Cauchy scale shared
by all chains. With correlated predictors, `m0` is a coordinate reference
effective count, not prediction degrees of freedom or a literal nonzero count.

`tpb_calibrate_phi(X, m0)` also works separately; its `scale_phi` can be passed
in `hyper_params` to `fullGibbs()`. It uses deterministic Sobol integration from
the existing GSL dependency and leaves the R random seed unchanged.
`integration_change` compares the full and half-sized integration grids at the
same scale; it is a diagnostic, not an error bound. Increase `n_points` (or
`calibration_points` in the pipeline) to check accuracy. No integration grids
are retained. Without `m0`, calibration is skipped and the supplied
`scale_phi` is used for the selected global prior.

`fullGibbs()` combines marginal MH with a beta-collapsed update of a, b and phi
in transformed local scales. This preserves the TPB model and hyperpriors while
improving mixing. Local variances and very small coefficient magnitudes are kept
on the log scale internally. Set `interweave = FALSE` to skip the additional
update; `tpb_full_pipeline()` accepts the same option through `...`.

Use `sampler = "conditional"` for the original full-conditional augmentation
(method 1). Method 3 uses `sampler = "marginal", interweave = TRUE`. Both support
Gamma or half-Cauchy hyperpriors for a and b. The default is now
`sqrt(phi) ~ half-Cauchy(0, scale_phi)`. phi is still the variance multiplier
in `beta_j | scales ~ N(0, phi * nu_j/lambda_j)`; the coefficient prior does
not include sigma squared. To reproduce the earlier global prior, set
`prior_type_phi = "hcauchy"` in `hyper_params`. For the default
`prior_type_phi = "sqrt_hcauchy"`, `scale_phi` is a standard deviation scale;
for `"hcauchy"` it is a variance scale. Old calibrated values must be recomputed.

The base marginal MH step jointly proposes log(a) and log(phi), then updates
log(b) separately with proposal standard deviation 0.1. This is the former
`bi_adaptive` update. Its two-dimensional proposal covariance is
adapted during burn-in and fixed afterward. The `proposal_type` and
`mh_step_a`, `mh_step_b`, `mh_step_phi` arguments have been removed; omit them
from calls written for earlier versions. The conditional sampler and the
additional interweaving block retain their slice updates.

`Gibbs_beta()` no longer accepts the unused `a` and `b` arguments. It receives
the prior variances directly, or forms them as `phi * nu / lambda`. The EB
sampler uses the `omega` argument for that multiplier; it does not pass the
marginal TPB phi without converting between the two parameterizations.

```r
hp <- list(prior_type_a = "hcauchy", prior_type_b = "hcauchy",
           scale_a = 1, scale_b = 1,
           prior_type_phi = "sqrt_hcauchy", scale_phi = 1)

fit1 <- fullGibbs(X, y, num_burnin = 5000, num_output = 5000,
                  sampler = "conditional", hyper_params = hp,
                  adapt_slice = TRUE)
fit3 <- fullGibbs(X, y, num_burnin = 5000, num_output = 5000,
                  hyper_params = hp, interweave_every = 3,
                  adapt_slice = TRUE)
```

`interweave_every = 3` runs the additional block every third iteration; the base
sampler runs every iteration. `adapt_slice = TRUE` adjusts the starting slice
width every 100 actual updates during burn-in, using twice the mean absolute
displacement. Widths are fixed after burn-in. The previous defaults remain
`interweave_every = 1` and `adapt_slice = FALSE`.

`tpb_full_pipeline()` uses iteration totals across chains. Four chains with
5,000 burn-in and 5,000 retained draws each require `num_iter = 40000`,
`num_warmup = 20000`, `num_chains = 4`. For half-Cauchy use `auto_find = FALSE`
and pass `hp`, with optional `m0` for calibration; `auto_find = TRUE` instead selects the Gamma hyperprior via the
existing reverse-logistic procedure. In sampler comparisons, run selection
once and pass the same resulting hyperparameters to both samplers.

UTPB keeps independent `eta_j ~ Uniform(0, 1)`, `a_j = eta_j`, `b_j = 1-eta_j`.
Its `tpb_fullyBayes_eta_uniform()` entrypoint accepts the same sampler,
interweaving interval and slice adaptation options. The additional UTPB block
uses joint MH for local logit(eta), a slice update for log(phi), and a beta
redraw. Use `store_locals = FALSE` when nu/lambda draws are not needed. The
conditional sampler's default still returns them; the marginal sampler does not.

Final slice widths are returned for reproducibility. The conditional sampler
does not compute marginal log densities; these fields are `NULL`, and its
MH acceptance rates are `NA` because its hyperparameters use slice updates.

Calibration and baseline calls

`tpb_calibrate_phi()` supports random Gamma or half-Cauchy shape parameters,
fixed shapes (`model = "fixed", hyper_params = list(a = ..., b = ...)`), and
UTPB (`model = "utpb"`). UTPB integrates the joint specification
`a = eta, b = 1-eta`, not two independent uniforms. The samplers for fixed
shapes and UTPB accept `m0` directly. Calibration runs once before sampling;
for several chains, calibrate once and pass `scale_phi` to each chain instead.
The mean target always averages over the chosen priors, including the random
global scale. The Gamma case still integrates a and nu analytically; the other
cases integrate nu analytically and average over the remaining shape prior.

For the four-way prior comparison, use Gamma or half-Cauchy hyperparameters
with `tpb_full_pipeline()`, each with and without `m0`. Keep the sampler and
noise model the same in all four calls. `m0 = NULL` uses the supplied
`scale_phi`; `m0` supplied solves for it from the prior mean effective count.

```r
# Strawderman-Berger: fixed TPB shapes, calibrated global scale.
sb <- tpb_fullyBayes_fixab(X, y, a = 1, b = 0.5,
                           num_burnin = 5000, num_samples = 5000,
                           m0 = m0, store_locals = FALSE)

utpb <- tpb_fullyBayes_eta_uniform(X, y, 5000, 5000, m0 = m0,
                                  sampler = "marginal", interweave_every = 3,
                                  adapt_slice = TRUE, store_locals = FALSE)

# Horseshoe with the PV reference rule, for standardized columns.
p <- ncol(X)
pv_scale <- m0 / (p - m0) / sqrt(nrow(X))
hs <- tpb_fullyBayes_fixab(X, y, a = 0.5, b = 0.5,
                           num_burnin = 5000, num_samples = 5000,
                           scale_phi = pv_scale, scale_with_sigma = TRUE,
                           store_locals = FALSE)

# IGG: original local prior, transferred to p regression coefficients.
igg <- tpb_fullyBayes_fixab(X, y, a = 1/p, b = 0.5 + 1/p,
                            num_burnin = 5000, num_samples = 5000,
                            fix_phi = TRUE, phi_val = 1, store_locals = FALSE)
```

For `scale_with_sigma = TRUE`, the fixed-shape sampler uses coefficient variance
`sigmaSq * phi * nu/lambda`, updates sigma squared accordingly, and returns
phi before multiplication by sigma squared. Thus its physical global standard
deviation is `sqrt(sigmaSq * phi)`, with conditional half-Cauchy scale
`sqrt(sigmaSq) * scale_phi`. Only this baseline option changes the noise scaling;
the main TPB and UTPB models retain their original variance convention.
The PV call deliberately leaves `m0 = NULL`: it uses the PV reference formula,
not our integrated mean calibration. If both `m0` and `scale_with_sigma = TRUE`
are used, the mean calibration cancels sigma squared from the shrinkage factor.

The IGG paper's inverse-Gamma shape is our b, and its Gamma shape is our a.
Its normal-means dimension is the number of coefficients, p in this regression
adaptation. The previous experiment driver had these shapes reversed and sampled
an extra global multiplier. The call above uses the corrected mapping and fixed
phi = 1. It retains our regression likelihood and noise prior; it is not a
literal reproduction of the paper's known-noise normal-means experiment.
See [Bai's dissertation, Eq. (2-5) and Section 2.4.1](https://ufdcimages.uflib.ufl.edu/UF/E0/05/24/20/00001/BAI_R.pdf)
and [Piironen and Vehtari, Eq. (16)](https://proceedings.mlr.press/v54/piironen17a/piironen17a.pdf).
