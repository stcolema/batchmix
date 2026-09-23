# batchmix 3.0.0

## Breaking (soft - old code still works, with a deprecation warning)

* **`fitBatchMix()` is now the package's main entry point**, replacing
  `runMCMCChains()` as the recommended way to fit a model: it runs several
  chains (fitting only one chain is not extra rigour, it is how
  non-convergence is detected at all), auto-tunes every Metropolis-Hastings
  proposal window by default, and reports rank-normalized split-Rhat/ESS and
  the best chain by BIC automatically. `runMCMCChains()` is now a thin
  deprecated alias for `fitBatchMix()` (forwards every argument unchanged,
  warns once per call) and will be removed in a future release.
* **The iteration-count argument has been renamed `R` -> `n_iter`** on every
  entry point (`runBatchMix()`, `batchSemiSupervisedMixtureModel()`,
  `fitBatchMix()`, `continueChain()`, `continueChains()`), matching common R
  MCMC convention (`nimble::runMCMC(niter = )`, `rjags`' `n.iter`) and
  avoiding the collision with `R`-as-in-"the R language" and with the
  correlation matrix `R` used elsewhere in this package's own `type =
  "MVN_LKJ"`/`"MVN_MIXED"` documentation. The old `R` argument still works
  this release (with a warning) when passed positionally in its original
  slot, or when every other argument in the same call is also named; mixing
  a named `R` with later positional arguments is not supported (there is no
  way to bind two names to the same argument slot - update such calls to use
  `n_iter` or name every argument).
* A message is now printed on `library(batchmix)` summarising the two
  changes above; silence it with
  `suppressPackageStartupMessages(library(batchmix))`.
* The dense, unordered parameter lists of `runBatchMix()`,
  `batchSemiSupervisedMixtureModel()` and `fitBatchMix()` have been
  regrouped (problem specification / MCMC control / proposal windows / prior
  hyperparameters / initial values / optional structural extensions), with
  matching section comments in the source. This is a pure reordering of
  already-named-only arguments - nothing that was ever called positionally
  moved position - but it does change `formals()`/autocomplete order.
* **`runBatchMix()` is now the single documented owner of every shared
  argument's description.** `batchSemiSupervisedMixtureModel()` and
  `fitBatchMix()` previously each carried their own copy (of varying
  completeness) of the same argument docs; they now inherit them from
  `runBatchMix()` (`@inheritParams`) instead of repeating them, so there is
  exactly one place to read about e.g. `mu_proposal_window`, and it cannot
  drift out of sync between the three. `batchSemiSupervisedMixtureModel()`'s
  docs now say plainly that it is the low-level engine underlying the other
  two (despite its name, it already handled the unsupervised case
  identically) and that most users should call `fitBatchMix()`/
  `runBatchMix()` instead.
* **Every Metropolis-Hastings proposal window, plus `auto_tune`/`n_burn`,
  is now bundled into a single `control = batchmixControl(...)` argument**
  on `runBatchMix()`, `fitBatchMix()` and
  `batchSemiSupervisedMixtureModel()` (Stan/`glmer`-style), replacing 11
  individual, rarely-touched arguments that previously sat among the main
  argument list. Prior/model-structure hyperparameters that are rarely
  touched for a different reason - because they change what is fitted, not
  how hard the sampler works to fit it (`rho`, `theta`, `m_scale`, `eta`,
  `a_gamma`/`b_gamma`, `batch_weight_prior` and its GP/partial-pooling
  hyperparameters) - deliberately stay as ordinary named arguments rather
  than joining `control`. The 11 individual arguments still work this
  release (with a warning) as deprecated aliases, merged into `control`;
  if both are supplied, `control` wins (also with a warning). See
  `?batchmixControl`.

## New features

* **Fitted models now print and summarise like a model object, not a raw
  list** (`runBatchMix()`/`batchSemiSupervisedMixtureModel()`/
  `continueChain()`/`processMCMCChain()` return an object of class
  `batchmix_fit`; `fitBatchMix()`/`continueChains()`/`processMCMCChains()`
  return a `batchmix_fit_list`). `print()` gives a short mclust/stanfit-style
  report (density type, dimensions, MCMC settings, mean acceptance rates and
  any attached convergence diagnostics) instead of dumping every sampled
  array to the console; `summary()` gives more detail (BIC, log-likelihoods,
  per-parameter-family acceptance rates). `assessConvergence()`'s return
  value (class `batchmix_convergence`) gets the same treatment, and its
  `format()` method is now the single place the Rhat/ESS/best-chain report
  text lives - `fitBatchMix()`'s automatic post-fit message uses it too, so
  the wording can no longer drift between the two. These are purely
  additive: every object is still a plain list/list-of-lists underneath, so
  existing `$`/`[[`-based code is unaffected.
* **Prior and posterior predictive checks**
  (`simulatePriorPredictive()`/`simulatePosteriorPredictive()`/
  `plotPredictiveCheck()`): the standard Bayesian-workflow bracket around a
  fit (BDA3 Ch. 6; Gelman et al. 2020, "Bayesian Workflow") - simulate
  datasets from the prior before fitting to sanity-check the (empirical-
  Bayes) hyperparameters, and from the fitted posterior afterwards to check
  the fit reproduces the observed data's features, including per-batch
  statistics specifically targeting whether the batch correction worked.
  Not yet supported for `type = "MVN_MIXED"` (documented as a scope
  boundary, not a silent gap).
* **Rank-normalized, folded, split-Rhat and bulk/tail ESS**
  (`rankNormalizedRhat()`), the diagnostic that superseded classical
  Gelman-Rubin R-hat (Vehtari, Gelman, Simpson, Carpenter & Burkner, 2021,
  *Bayesian Analysis* 16(2)), computed automatically by `fitBatchMix()`
  whenever `n_chains >= 2`, on the permutation-invariant complete-data
  log-likelihood trace (raw cluster parameters are not label-switching
  invariant). `assessConvergence()`/`getBestChain()` expose this
  programmatically; the best chain (by post-burn-in mean BIC) is attached
  as `attr(., "best_chain")` and propagated through
  `processMCMCChains()`/`continueChains()` non-breakingly (via attributes,
  not a return-shape change).
* **`generateInitialLabels()` seeds unsupervised fits from a k-means
  partition of the data by default** (when `X` is available and no items
  are fixed), instead of a random draw from the stick-breaking prior. A
  random draw is completely uninformed by the data; for well-separated
  clusters this made it easy for a chain to start in (and then, since
  single-site Metropolis-Hastings essentially never crosses a low-density
  valley in a realistic run length, get permanently stuck in) a poor local
  optimum - diagnosed as the cause of a reported Rhat that got *worse*
  rather than better with a longer chain (different chains stuck in
  different, persistently-worse-fitting optima, not label-switching: no
  column-swap was observed in any chain's cluster-mean trace, and the
  auto-tuned proposal windows were unremarkable). This does not fully
  eliminate the risk (k-means itself is not guaranteed to find the best
  partition) - `fitBatchMix()`'s Rhat warning now says so explicitly and
  suggests more chains rather than only a longer run.
* **`plotBatchCorrection(..., show_raw = FALSE)`**: the raw-vs-corrected
  comparison plot can now show just the corrected panel.
* **Missing entries in `X` (`NA`/`NaN`) are now supported by every sampler
  type**, not only `type = "MVN_MIXED"`: `MVN`, `MVT` and `MVN_LKJ` each
  gained an `updateLatentData()` Gibbs step generalising the construction
  `mvnSamplerMixed::updateLatentData()` already used (Albert & Chib, 1993;
  Dunson, 2000) - every missing entry of a row is redrawn from its full
  conditional given the row's other (observed or already-imputed) entries
  and the current cluster/batch parameters, via the precision matrix
  (`cov_comb_inv`) each sampler already maintains, once per MCMC sweep,
  before that sweep's allocation/parameter updates - never imputed once and
  held fixed. `MVT` additionally draws a transient per-item Gamma weight
  each sweep (the standard Gaussian-scale-mixture representation of the
  multivariate t; Liu & Rubin, 1995) so the imputation targets the correct
  Student-t - not Gaussian - full conditional; the weight is not retained
  state, so every other MVT update still targets the true marginal
  t-density on the completed row, exactly as when `X` has no missing data.
  This is correct under MCAR, and under MAR given the modelled cluster/
  batch structure and the other observed columns - not MNAR (e.g. values
  missing because they are unusually large/small), which no sampler here
  attempts to model. Every constructor's empirical-Bayes prior setup
  (previously raw `mean(X)`/`cov(X)`, undefined under `NA`) now uses a
  column-mean-imputed copy of `X` for that one-off calculation only (see
  `imputeColumnMeans()` in `src/genericFunctions.cpp`), matching what
  `mvnSamplerMixed` already did. The completed data is tracked separately
  from the user's own `X` (never modified) as a proper per-iteration trace,
  `mcmc_output$latent_data`, combined correctly across
  `continueChain()`/`continueChains()`. `simulatePriorPredictive()` also no
  longer rejects a missing-data `X` (it only ever needed `X`'s column
  mean/covariance, which are now computed the same NaN-tolerant way).
  `type = "MVN_MIXED"` remains the only type additionally modelling
  censored and binary/probit columns.

## Internal

* Collapsed the diamond/virtual-multiple-inheritance C++ class hierarchy
  (`sampler`/`semisupervisedSampler` and the `*Predictive` wrapper classes)
  into a single inheritance chain per sampler family, verified
  behaviour-preserving via bit-exact golden-master tests. Deduplicated the
  four near-identical MCMC driver loops into one shared templated function
  (`runSemisupervisedSampler()` in `src/runSampler.h`).
* Fixed systemic roxygen2 documentation bugs found while auditing every
  vignette against the current API: 108 tag lines in `genericFunctions.cpp`
  missing the `@` prefix roxygen requires, a multi-name `@param` list in
  `runBatchMix()` silently truncated by a line wrap, and doc-comments
  written on the (dead, from roxygen2's perspective) `.h` declaration
  instead of the `.cpp` `// [[Rcpp::export]]` definition.



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

* Fixed `fitBatchMix()` rejecting every argument forwarded through its own
  documented `...` (e.g. `include_interaction = TRUE`,
  `r_proposal_window = 0.05`) with an "unused argument" error - the
  internal check backing the `R` -> `n_iter` shim was written for
  `runBatchMix()`/`batchSemiSupervisedMixtureModel()`, whose `...` accepts
  only the deprecated `R`, and was reused as-is for `fitBatchMix()`, whose
  `...` has a much wider, documented contract.
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
* **Four vignettes**, each a self-contained worked example following the
  same loop (simulate against a known ground truth, fit several chains,
  check MCMC diagnostics via `print()`/`summary()`, and verify recovery),
  replacing what were six more narrowly-scoped and partly-overlapping
  documents:
    - `vignettes/batchmix_workflow.Rmd` ("Getting started"): the full
      loop end to end on `type = "MVN"`, including the prior/posterior
      predictive check bracket (`simulatePriorPredictive()`/
      `simulatePosteriorPredictive()`).
    - `vignettes/covariance_models.Rmd`: comparing `"MVN"`, `"MVT"` and
      `"MVN_LKJ"` - cluster-count comparison via BIC, MVT's robustness to
      heavy-tailed data, and LKJ correlation recovery.
    - `vignettes/batch_weight_priors.Rmd`: `"global"`, `"partial_pooling"`
      and `"gp"` batch-weight priors on a simulated six-batch scenario with
      a smooth assay drift and a genuine step change in per-batch class
      proportion.
    - `vignettes/probit_missing_censored.Rmd`: `type = "MVN_MIXED"` for
      binary/probit columns and censored entries.
* Added `tests/testthat/test-new-features.R`, regression tests for each of
  the bugs above and for the new features' basic recovery properties.
