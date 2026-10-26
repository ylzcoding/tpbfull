EM-within-Gibbs for re-parameterized three-parameter beta model

`fullGibbs()` combines marginal MH with a beta-collapsed update of a, b and phi
in transformed local scales. This preserves the TPB model and hyperpriors while
improving mixing. Local variances and very small coefficient magnitudes are kept
on the log scale internally. Set `interweave = FALSE` to skip the additional
update; `tpb_full_pipeline()` accepts the same option through `...`.

Use `sampler = "conditional"` for the original full-conditional augmentation
(method 1). Method 3 uses `sampler = "marginal", interweave = TRUE`. Both support
Gamma or half-Cauchy hyperpriors for a and b; phi itself has a half-Cauchy prior.

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
           scale_a = 1, scale_b = 1, scale_phi = 1)

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
and pass `hp`; `auto_find = TRUE` instead selects the Gamma hyperprior via the
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
