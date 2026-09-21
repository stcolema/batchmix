# batchmix 2.1.0

## New features

* **Batch x cluster interaction term** in the mean (`include_interaction`
  argument to `batchSemiSupervisedMixtureModel()`/`runBatchMix()`): a
  partial-pooling `gamma_{k,b} ~ N(0, tau2_interaction)` term, confined to
  the sum-to-zero subspace (every row and column of `gamma[p, , ]` sums to
  zero) so it is identifiable against the existing cluster mean/batch shift
  main effects rather than confounded with them (the classical two-way
  ANOVA main-effect/interaction problem). Available for every sampler type.
* **Batch-specific mixture weights** (`batch_weight_prior` argument, one of
  `"global"`, `"partial_pooling"`, `"gp"`): each batch can now have its own
  mixture weight vector instead of every batch sharing one. `"global"`
  reproduces the original behaviour exactly. `"partial_pooling"` gives each
  batch an independent, exchangeable additive-log-ratio weight shrunk
  towards a common conjugate Gibbs-sampled mean/variance
  (`eta_{b,j} ~ N(mu_j, tau2_j)`) - for batches whose proportions differ but
  have no known order or spatial/temporal structure. `"gp"` instead links
  the free ALR coordinates across batches with a Gaussian process prior
  over an ordered/spatial `batch_coordinates` covariate - for batches
  collected over time or space whose cluster proportions are expected to
  vary smoothly. GP hyperparameters (`gp_tau2`, `gp_length_scale`) can be
  fixed or sampled (`sample_gp_hyperparameters`).
* **Auto-tuning of Metropolis-Hastings proposal windows**
  (`auto_tune`/`n_burn` arguments, `auto_tune = TRUE` by default): every
  proposal window is adapted via Robbins-Monro diminishing adaptation
  during burn-in and then frozen, rather than requiring hand-tuned windows.
  Final tuned windows are returned in the output (`final_*_proposal_window`
  fields) and reused by `continueChain()`/`continueChains()`.
* The LKJ/separation-strategy sampler (`type = "MVN_LKJ"`: LKJ prior on the
  cluster correlation matrix, log-normal marginal scales, decoupling
  correlation and scale beliefs) and the mixed continuous/binary/missing/
  censored-data sampler (`type = "MVN_MIXED"`) are now available from the
  main `batchSemiSupervisedMixtureModel()`/`runBatchMix()` API, including
  the semi-supervised (`fixed`) path, rather than only via low-level direct
  constructor calls.

## Bug fixes

* Fixed a Metropolis-Hastings sign error in the batch scale
  (`batchScaleMetropolis()`) and marginal-scale (`sigmaMHStep()`) updates
  across every sampler variant: the forward and reverse proposal densities
  for the asymmetric Gamma random-walk proposal were assigned to the wrong
  side of the acceptance ratio, silently biasing the sampled batch
  scale/marginal-scale values (verified by simulation against a known
  Gamma(5, 2) target: the bug recovered a posterior mean of ~1.4 instead of
  the true 2.5). `mvtSampler::clusterDFMetropolis()` already had the
  correct orientation and was used as the reference to fix the rest.
* Fixed numerically catastrophic conditioning of the GP covariance matrix
  underlying `correlated_weights` for realistic length-scale choices (a
  fixed jitter of 1e-6 gave a condition number of ~1.7e7 in one verified
  case, corrupting every downstream calculation using its inverse); the
  jitter now scales with the GP marginal variance.
* Fixed `clusterMeanMetropolis()`/`batchShiftMetorpolis()` silently
  dropping the interaction term's contribution to `mean_sum` whenever the
  cluster mean or batch shift was updated with `include_interaction = TRUE`.
* Fixed `updateBatchCorrectedData()` only ever subtracting the batch shift
  `m_b` and never the interaction term `gamma_{k,b}`, so "batch-corrected"
  data for any (cluster, batch) cell with a real interaction effect still
  carried a leftover batch-specific offset the size of that interaction.
* Fixed `continueChain()`/`continueChains()` silently reverting a
  continued chain to `batch_weight_prior = "global"` and
  `include_interaction = FALSE` regardless of what the original chain
  used, and never forwarding the associated hyperparameters
  (`batch_coordinates`, `a_gamma`/`b_gamma`, `pp_tau2_shape`/`pp_tau2_rate`/
  `pp_mu_prior_sd`, etc.) - a continuation was silently changing what model
  was being fitted, not just resuming it. Also fixed the new weight-prior
  and interaction trace fields (`w_batch`, `eta_alr`, `gp_tau2`,
  `gp_length_scale`, `pp_mu`, `pp_tau2`, `gamma`, and their acceptance
  rates) not being combined with the original chain's samples when
  `keep_old_samples = TRUE`, so they were silently truncated back down to
  only the newly-sampled segment - including for the default
  `batch_weight_prior = "global"` case, where `w_batch`/`eta_alr` are still
  returned (just constant across batches) and so still need combining.

## Other changes

* Removed `src/RESHUFFLE/`, an earlier, uncommitted, non-compiling
  architecture experiment superseded by the above.
* Extended `vignettes/batchmix_workflow.Rmd` into a full worked example:
  simulating data with a visible batch effect, visualising it with and
  without the batch effect, fitting the model, checking convergence
  (acceptance rates, likelihood traces, and a Gelman-Rubin R-hat check),
  and checking parameter recovery and the batch correction itself against
  the simulation's ground truth. Added lighter data-visualisation sections
  to the other three vignettes.
* Added `vignettes/batch_weight_priors.Rmd`, comparing `"global"`,
  `"partial_pooling"` and `"gp"` on a simulated six-batch ELISA-style
  scenario with a smooth assay drift and a genuine step change (an
  "outbreak") in the true per-batch positive-class proportion, collected at
  irregular time points: visualising the data with and without the batch
  effect, fitting all three via `runMCMCChains()`, MCMC diagnostics
  (acceptance rates, likelihood traces, Gelman-Rubin R-hat, `continueChains()`),
  recovered vs true per-batch proportions, batch-corrected data by
  predicted class against ground truth, and classification accuracy by
  batch. Semi-supervised labels are chosen as the most extreme observed
  values per batch and class, not a random subset.
* Added `tests/testthat/test-new-features.R`, regression tests for each of
  the bugs above and for the new features' basic recovery properties.
