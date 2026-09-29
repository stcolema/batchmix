# batchmix 3.0.0

## Breaking changes to defaults and results (read first)

* **Partial pooling is now the default.** `batch_weight_prior` defaults to
  `"partial_pooling"` (was `"global"`; a single batch still resolves to
  `"global"`), and `sample_s_scale` defaults to `TRUE` (was `FALSE`), so the
  concentration of the batch-scale prior is estimated instead of fixed.
  `sample_s_scale`, `s_scale_shape` and `s_scale_rate` are now also exposed
  on `runBatchMix()` and `fitBatchMix()` (previously only on
  `batchSemiSupervisedMixtureModel()`). Restore the earlier model with
  `batch_weight_prior = "global", sample_s_scale = FALSE`. Consequences: (i)
  `alpha`/`concentration` are unused under partial pooling, but surplus
  components are still largely emptied because the partial-pooling weight
  prior is diffuse (`pp_mu_prior_sd = 10` on the logit scale: pairwise logit
  differences have standard deviation about 14, that of a Dirichlet with
  concentration about 0.1): with 3 true clusters and `K_max = 8`, about 3.8
  components stayed occupied, against 4.1 under the old global
  Dirichlet(1/`K_max`) default and 7.3 under Dirichlet(1). Sharpening
  `pp_mu_prior_sd` to 2.5 gave 4.8 and to 1 gave 7.9 (see
  `inst/experiments/overfitted_weight_priors.R`). This is empirical, not a
  theorem, and from one setting. (ii) The same diffuseness makes weights
  close to degenerate a priori, determined by the data. (iii) **The weight
  prior is now exchangeable over clusters.** Batch weights are the softmax of
  `K` logits with the same prior for every class (population mean
  constrained to sum to zero, class variances `tau2_k`), instead of `K - 1`
  log-ratios against the last class as reference; the reference
  parameterisation made the last cluster the dominant one with prior
  probability 0.008 (against 0.14 for each other cluster at `K_max = 8`). The
  common per-batch shift that softmax cannot see is updated by its own exact
  Gibbs step. This applies to `"gp"` too (`K` independent GPs with sum-to-zero
  intercepts). Outputs change shape and name accordingly: `eta_alr` is now
  `eta_logit` (`B x K x iterations`), and `pp_mu`, `pp_tau2` and `gp_beta` have
  `K` rows. Every fit's fixed-seed output changes.
* `simulatePriorPredictive()` follows the default model: batch-specific
  weights under partial pooling (`w_batch` in `params`), `rho` drawn from its
  hyperprior when `sample_s_scale = TRUE`, and `m_scale = NULL` (the default)
  draws the batch-shift scale from its InvGamma(3, 1) hyperprior instead of
  fixing it at 0.01, matching the fitting functions.
* The startup message now lists these changes and points to the vignettes,
  which were updated throughout (see the items below).

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
* **`simulatePriorPredictive()`/`simulatePosteriorPredictive()` now support
  `type = "MVN_MIXED"`.** Every column simulates the same shared latent
  Gaussian draw as `"MVN_LKJ"`; binary/probit-linked columns
  (`column_type == 1`) are then thresholded at 0 (Albert & Chib, 1993),
  matching `sigma`/`S` being fixed at 1 for those columns in the real
  sampler. Every replicate is fully observed (no `NA`/censoring) - the
  model's own assumption is that missing/censored cells are draws from
  exactly the same distribution as every other cell, so simulating them
  like any other cell is the model-consistent replicate.
  `plotPredictiveCheck()` gained a `censor_code` argument that instead
  excludes the real data's censored cells from the comparison (a recorded
  censoring bound is not the true value to compare a free replicate draw
  against), and now warns if asked for a `style = "density"` check on an
  apparently-binary column (`style = "statistic"` with `statistic = mean`,
  the proportion of 1s, is the appropriate check there).
* **`predictNewBatch()`**: posterior-predictive quantities for a batch not
  seen during fitting, for every model `type` and every `batch_weight_prior`.
  With no new data, draws the new batch's own class weight and location
  shift/scale straight from their fitted priors - under `batch_weight_prior
  = "gp"` this is a genuine Gaussian process prediction at a new
  `new_batch_coordinate` (e.g. a future collection time), and under
  `"partial_pooling"` a draw from the estimated population distribution of
  batch weights (the standard "predict a new group" logic for a
  hierarchical model). With the new batch's own data supplied (`X_new`),
  additionally classifies it by composition sampling: for each retained
  posterior draw, the new batch's shift/scale/weight and item allocations
  are drawn from their conditional posterior given that draw's (fixed)
  cluster parameters and `X_new`'s own likelihood, reusing the training
  sampler's own Metropolis/Gibbs updates restricted to the one new batch -
  a proper posterior predictive by composition (Rubin, 1987), not a
  plug-in point estimate. See `?predictNewBatch` for the full statistical
  design.
* **Partial pooling for the batch-scale prior's concentration**
  (`sample_s_scale` argument to `batchSemiSupervisedMixtureModel()`/
  `runBatchMix()`). The batch shift prior already had a genuine
  partial-pooling hyperparameter (`lambda_2`, via `sample_m_scale`); the
  batch scale prior `S_b ~ S_loc + InvGamma(rho, theta)` previously fixed
  both `rho` and `theta` as user-supplied constants, never re-estimated
  from the batches actually fitted. An additive/multiplicative shift-scale
  decomposition is only identifiable up to one free hyperparameter per
  batch effect (a uniform rescaling of every batch's `S_b` is exactly
  compensated by the inverse rescaling of the cluster covariance), so only
  the InvGamma's concentration `rho` - never its mean - can be freely
  estimated; `theta` is kept in lock-step (`theta = s_scale_prior_mean *
  (rho - 1)`) so the prior's mean stays fixed at its original,
  user-supplied value while `rho` alone adapts to how similar or different
  the fitted batches' own scales turn out to be. Off by default
  (`sample_s_scale = FALSE` reproduces the previous fixed-prior behaviour
  exactly); when `TRUE`, `rho` is updated by a symmetric log-random-walk
  Metropolis step (`s_scale_shape`/`s_scale_rate` set its own Gamma prior,
  `s_scale_proposal_window` its auto-tuned proposal window - see
  `?batchmixControl`). `predictNewBatch()`'s composition sampling for a new
  batch's scale now conditions on each posterior draw's own `rho` when the
  original fit used `sample_s_scale = TRUE`.

## Bug fixes

* Fixed every `sampleSemisupervised*()` C++ driver building its return value
  with a single `Rcpp::List::create(Named(...) = ..., ...)` call of 39-42
  arguments - `List::create()`'s generated overloads only go up to 20
  arguments in Rcpp versions before its May 2024 move to a genuine C++11
  variadic template (first released as Rcpp 1.0.13), so this failed to
  compile against any such Rcpp - including the package's own declared
  minimum, `Rcpp (>= 1.0.5)`. Replaced with incremental `out["key"] <-
  value` assignment, which has no such limit on any Rcpp version.
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
* Fixed auto-tuning moving the cluster-covariance (Wishart df), batch-scale,
  marginal-scale (`sigma`) and t degrees-of-freedom proposal windows in the
  wrong direction. The R layer hands C++ the reciprocal of these windows
  (larger = tighter proposal), but the Robbins-Monro update treated them as
  step sizes, so high acceptance widened the acceptance-driving window
  instead of loosening the proposal. They now use
  `robbinsMonroUpdateReciprocal()`; the covariance df is also floored at
  `P`, below which the Wishart proposal is invalid. Acceptance rates feeding
  adaptation are now averaged over occupied clusters only.
* Fixed `MVN_LKJ`/`MVN_MIXED` marginal-scale (`sigma`) prior: the Gamma
  random walk on `sigma` was scored against the log-normal density without
  its `1 / sigma` Jacobian, so the chain targeted `LogNormal(beta + xi^2,
  xi)` (median ~2.3 at the defaults) instead of `LogNormal(beta, xi)`
  (median ~0.85).
* Fixed `r_log_det` being left at zero until the first accepted correlation
  proposal, so the first `rMHStep()` decisions scored the LKJ prior term
  `(eta - 1) * log det(R)` incorrectly for every `eta != 1`.
* Fixed empty-cluster forced prior draws being counted as Metropolis
  acceptances (`cov`/`mu`/`r`/`sigma`/`t_df` counters), which inflated both
  the reported `*_acceptance_rate` values and the adaptation signal. Also
  fixed the `MVN_LKJ`/`MVN_MIXED` empty-cluster draw order: `mu` is now drawn
  after `R`/`sigma` so that `(mu, Sigma)` follows the joint prior.
* **Changed what `batch_corrected_data`/`inferred_dataset` are.** They were
  `(x - mu_k - m_b) / sqrt(S_b) + mu_k`, which is not the inverse of the
  fitted batch model: the model inflates only the diagonal of the cluster
  covariance (`Sigma_kb = Sigma_k + diag((S_b - 1) * diag(Sigma_k))`), so
  dividing the residual by `sqrt(S_b)` on both sides restored the marginal
  variances but shrank between-feature correlations by `1 / sqrt(S_p S_q)`.
  They are now the posterior mean of each item's batch-free signal,
  `mu_k + Sigma_k Sigma_kb^-1 (x - mu_k - m_b)` (with the interaction term
  subtracted along with `m_b` when in use; exact for MVN, and for MVT
  because the latent scale cancels), computed by the new internal
  `batchCorrectedPosteriorMean()`. This is a denoising estimate: within-
  cluster deviations are shrunk towards `mu_k`, so the corrected data are
  less variable than the true batch-free data. Also corrected the
  `probit_missing_censored` vignette, which stated the congruence form
  `Delta_b Sigma Delta_b` rather than the diagonal inflation the code uses.
* **`salso` is no longer a dependency.** It needed a Rust toolchain to build
  from source, which is a barrier on clusters and managed laptops without
  admin rights. `predictFromMultipleChains()` now uses `minVI(method =
  "greedy")` on the posterior similarity matrix: the best of the
  `avg`/`comp`/`draws` candidates refined by greedy single-item relocation
  under the same Jensen lower bound to the expected Variation of Information
  (new `minVIGreedyRefine()`, in C++). On simulated planted-partition
  problems it matched `salso`'s exact expected VI. The trade-off: `minVI()`
  works from the `N x N` similarity matrix (`8 N^2` bytes), so it is not
  suitable for tens of thousands of items, where `salso`'s PSM-free search
  would be. `nCores` is now ignored. The covariance-models vignette now uses
  `minVI()` (VI loss) in place of `salso(loss = "binder")`. Also fixed
  `VI.lb()` rejecting a clustering that carries attributes (as `minVI()`'s
  output does).
* **New `calcBICM()`, and BICM reported by `summary()` and
  `assessConvergence()` (`chain_bicm`).** The stored per-iteration `BIC`
  plugs the log-likelihood at a single posterior draw into the BIC formula;
  a draw sits below the maximum-likelihood fit, so this over-penalises
  complexity and is noisy. BICM (Raftery et al., 2007) uses the best
  post-burn-in draw instead. The per-iteration `BIC` is kept as a trace, and
  chain selection is unchanged. BICM shares the missing-data caveat
  documented at `?calcBICM`: the working (imputed/latent) likelihood is not
  the observed-data likelihood.
* Efficiency: the cluster-mean random walk is now shaped by each cluster's
  covariance (`mu' = mu + w chol(Sigma_k) z`; `mu_proposal_window` is now in
  units of the cluster's spread); the partial-pooling batch-weight update
  proposes each batch's entry separately rather than as one B-dimensional
  block (`eta_acceptance_rate` is now a per-move rate); and the interaction
  update scores the change in log-likelihood incrementally instead of
  re-evaluating every item's full quadratic form twice per feature. The
  GP-correlated weight update keeps its block proposal, since its prior
  correlates the batches. Each change is checked against the exact target
  distribution in `test-mcmc-kernel-fixes.R`.
* Robustness: a label draw can no longer index out of range when the
  cumulative allocation probabilities round below 1; `sample_s_scale = TRUE`
  with `rho <= 2` now errors instead of silently leaving `rho` frozen at its
  initial value; dead empirical-Bayes code and the constructors' console
  output were removed.
* Numeric output changes for a fixed seed: every golden fixture whose run has
  an empty cluster or uses `MVN_LKJ`/`MVN_MIXED` is regenerated (see
  `tests/testthat/test-mcmc-kernel-fixes.R` for the prior-recovery tests).
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
* Fixed `sampleMScalePosterior()` (the conjugate Gibbs update for the
  batch-shift variance hyperparameter `lambda_2`) computing the InvGamma
  posterior rate as `sum(m^2) / (4 * delta_2)` instead of the correct
  `sum(m^2) / (2 * delta_2)` - an erroneous extra factor of 0.5, live every
  sweep by default (`sample_m_scale = TRUE`), that systematically shrunk
  `lambda_2` and so over-shrunk every batch shift `m_b` toward zero more
  than the model's own prior warranted. Present identically in
  `mvnSamplerSeparationStrategy` (and hence `mvnSamplerMixed`, which
  inherits it).
* Fixed `sampleMPrior()` drawing the prior for the batch shift `m` as
  `N(mean, precision)` instead of `N(mean, 1/precision)`
  (`batch_shift_prior_precision` is a precision, so its standard deviation
  is the inverse square root, not the precision itself) - with typical
  hyperparameters this drew the initial `m` roughly three orders of
  magnitude too diffuse. Only affects `sampleFromPriors()`'s one-off initial
  draw (the ongoing Metropolis-Hastings step already used the correct
  density), so burn-in started from a badly-scaled point rather than the
  stationary posterior itself being affected. Present identically in
  `mvnSamplerSeparationStrategy`.
* Fixed `mvnSamplerSeparationStrategy`/`mvnSamplerMixed`'s `sigmaMHStep()`
  for a persistently empty cluster (`N_k(k) == 0`): despite the comment
  ("sample from the prior distribution"), it actually proposed an
  uncorrected, force-accepted random walk off the cluster's current
  (possibly stale) marginal scale `sigma`, rather than a fresh draw from the
  LogNormal(beta, xi) prior - an unbounded drift for any component with no
  members, inconsistent with the (correct) empty-cluster handling already
  used for `mu`/covariance elsewhere in the same file.
* Fixed the semi-supervised observed-data log-likelihood (and hence BIC)
  marginalising over every component even for items with a known (`fixed`)
  label, instead of using the joint density at the known label directly -
  this systematically inflated `observed_likelihood`/`BIC` for any
  semi-supervised fit, biasing `getBestChain()`'s BIC-based chain selection
  and any BIC comparison across models with different `fixed` vectors.
* Fixed `calcBIC()` (every sampler variant) not counting the extra free
  parameters contributed by the opt-in cluster x batch interaction term
  (`include_interaction = TRUE`) or a batch-specific weight prior
  (`batch_weight_prior = "partial_pooling"`/`"gp"`), biasing BIC comparisons
  in favour of the richer model whenever either feature is used (both
  contribute exactly 0 extra parameters in the default configuration, so
  ordinary fits are unaffected).
* Fixed `processMCMCChain()` computing every cluster-indexed point estimate
  (`mean_est`, `cov_est`, `mean_sum_est`, `cov_comb_est`, `weights`,
  `t_df_est`, `gamma_est`, `w_batch_est`, `allocation_probability`) by
  averaging/taking the median of the raw per-iteration arrays directly by
  component index, with no correction for label switching - only valid when
  every component happens to be permanently anchored by a fixed
  (semi-supervised) item, silently wrong otherwise (any unsupervised fit, or
  any semi-supervised fit with `K_max` above the number of labelled
  classes). Added `relabelChain()` (Equivalence Classes Representatives via
  the Hungarian algorithm on each iteration's overlap with a reference
  labelling, `clue::solve_LSAP`) and call it, after burn-in and before any
  point estimate, from both `processMCMCChain()` and `calcAllocProb()`.
* Fixed `processMCMCChain()` withholding `allocation_probability`/`prob`/
  `pred` entirely for a fully unsupervised fit (only ever computing them for
  a semi-supervised one) - once `relabelChain()` (above) makes an
  unsupervised fit's relabelled `alloc` just as well-defined as a
  semi-supervised one's, there is no reason a point-estimate predicted
  label shouldn't be available for an unsupervised fit too (e.g. for
  comparing recovered clusters against a known ground truth on simulated
  data, as `vignette("covariance_models")` does).
* Fixed `processMCMCChain()`'s `mean_sum_est`/`cov_comb_est` extracting the
  wrong (cluster, batch) column block whenever `K_max > 1` and `B > 1`:
  `mean_sum`/`cov_comb` are stored with the cluster index varying slowest
  (column/block position `(k - 1) * B + b`), but the extraction pulled
  contiguous blocks as if the batch index varied slowest, silently mixing up
  which cluster's estimate got paired with which batch's.
* Fixed `processMCMCChain()`/`calcAllocProb()` applying no burn-in at all
  (dropping every sample instead of none) whenever the effective burn-in
  was 0 (e.g. any `0 < burn < thin`): the fix for the earlier `seq(1, 0) ==
  c(1, 0)` bug (see above in this same pass) introduced `-seq_len(0) ==
  -integer(0)`, and indexing with an empty vector in R - negated or not -
  selects nothing, not everything (`(1:5)[-integer(0)]` is `integer(0)`, not
  `1:5`); every burn-in trim in both functions now guards this case
  explicitly.
* Fixed `processMCMCChain()`'s `batch_corrected_data`/`alloc` trimming
  silently collapsing to one fewer dimension than every other trimmed
  quantity (breaking `point_estimate_method = "mean"` outright, and
  returning the wrong shape/un-averaged values for `"median"`) whenever
  `P == 1` or `K == 1` respectively, from a missing `drop = FALSE`.
* Fixed `minVI(..., max.k = <value>)` throwing `object 'k_inds' not found`
  for every method except `"draws"` - `k_inds` was only ever assigned inside
  the `is.null(max.k)` branch.
* Fixed `prepareInitialParameters()` requiring *both* dimensions of an
  initial means/batch-shift/batch-scale matrix to be wrong before rejecting
  it (`&` where `|` was meant), so e.g. a matrix with the right number of
  columns but wrong number of rows silently passed validation.
* Fixed `generateBatchData()`/`generateBatchDataMVT()`/
  `generateBatchDataVaryingRepresentation()` (and the hand-written simulator
  in `vignette("batch_weight_priors")`) applying `batch_scale` as a
  standard-deviation multiplier (`Var = sd^2 * batch_scale^2`) when the
  fitted model applies it linearly to the variance (`Var = sd^2 *
  batch_scale`) - simulating with a given `batch_scale` and fitting the
  model back to it recovered a `scale_est` roughly the square of the value
  used to generate the data, not the value itself.
* Fixed `simulatePriorPredictive()`'s batch-shift draw using the stale,
  pre-fix `1 / (delta_2 * lambda_2)` scale instead of the corrected
  `sqrt(delta_2 * lambda_2)` (see the `sampleMPrior()` fix above) - this R
  mirror of the C++ prior was not updated when the C++ itself was fixed,
  so the prior predictive check's simulated batch shifts were roughly two
  orders of magnitude too diffuse.
* Fixed `simulatePriorPredictive()`/`simulatePosteriorPredictive()`
  erroring or silently misbehaving for two ordinary configurations:
  `diag(cov_comb) <- ...` failed outright whenever `P == 1` (a univariate
  dataset - `cov[, , k]` collapses to a bare scalar, not a 1x1 matrix, for
  which `diag<-` doesn't work), and `simulatePosteriorPredictive()`'s
  `mean_sum`/`cov_comb` extraction collapsed to the wrong shape whenever
  `K_max * B == 1` (a single cluster and a single batch), from the same
  class of missing-`drop = FALSE` indexing gotcha fixed elsewhere in this
  pass.
* Fixed a crash reported against `type = "MVN_LKJ"` (and, by inheritance,
  `"MVN_MIXED"`): `inv_sympd(): matrix is singular or not positive
  definite`, and, once that was fixed, `Mat::operator(): index out of
  bounds`. Root-caused via a gdb backtrace on a real crashing run: a
  random-walk proposal for the correlation matrix `R`, reparameterised via
  a Cholesky/partial-correlation transform, is guaranteed PD in exact
  arithmetic for any candidate - but not kept away from the boundary of the
  PD cone, and a proposal legitimately close to that boundary can be PD in
  exact arithmetic yet numerically singular in floating point (confirmed:
  eigenvalues down to ~1.4e-4, smallest nominally -1.3e-16 after
  symmetrising - a real, if rare, floating-point edge case, not a coding
  error in the reparameterisation itself). Two compounding bugs followed:
  (1) every `inv_sympd(proposed_cov)`-style call across
  `rMHStep()`/`sigmaMHStep()`/`batchScaleMetropolis()`/`sampleCovPrior()` in
  both `mvnSamplerSeparationStrategy.cpp` and `mvnSamplerMixed.cpp` used the
  throwing form uncaught, crashing on a numerically-degenerate proposal
  instead of simply rejecting it (the textbook-correct treatment of a
  proposal whose target density is numerically undefined); (2)
  `rMHStep()`'s regular branch called the non-throwing, output-parameter
  form of `arma::chol()` but never checked its boolean return value - on
  failure that form leaves the output matrix empty (0x0) rather than P x P,
  and the very next line indexed it assuming P x P. Every such call site is
  now guarded: a numerically-degenerate Metropolis-Hastings proposal is
  treated as an automatic reject (matching the `next`/`continue` pattern
  already used elsewhere in this file for degenerate individual-parameter
  proposals), and a numerically-degenerate *prior* draw (in
  `sampleCovPrior()`, and the empty-cluster branches of `rMHStep()`/
  `sigmaMHStep()`) is redrawn (a negligible-measure rejection of the
  pathological tail, not a bias on the prior) rather than either crashing
  or silently accepting a broken matrix.
* Fixed `relabelChain()` (and hence `processMCMCChain()`/`calcAllocProb()`)
  crashing on real (as opposed to this function's own synthetic test
  fixtures) MCMC output, via `R CMD check --run-donttest` on
  `plotBatchCorrection()`'s own example. Two compounding bugs: (1)
  `mcmc_output$samples` is 0-indexed (0..K_max-1) - the raw C++ sampler's
  own convention - but `relabelChain()` was written and tested assuming
  1-indexed labels, so `perm[samples[t, ]]` silently dropped any item
  whose label was 0 (R indexing with 0 selects nothing, rather than
  erroring), corrupting `samples`' length; (2) `mcmc_output$t_df`, read via
  `$` (which does partial matching), silently matched
  `t_df_proposal_window` for a non-MVT fit instead of returning `NULL` -
  there is no field literally named `t_df` on an MVN/MVN_LKJ/MVN_MIXED fit
  - so `has_t_df` was wrongly `TRUE` and every downstream `t_df` access
  then operated on a bare scalar instead of an (n_saved x K) matrix. Fixed
  by working in 1-indexed space internally (converting at the function's
  boundary) and switching every `mcmc_output$field` read in this function
  to `mcmc_output[["field"]]` (exact match only, returns `NULL` rather than
  guessing when a field is genuinely absent) respectively. The three
  synthetic test fixtures for this function were also corrected to use
  0-indexed labels, matching real sampler output, since the original
  (1-indexed) fixtures were self-consistent but not representative and is
  exactly what let both bugs through in the first place.

## Other changes

* `R CMD check` is now clean (0 errors, 0 warnings, 0 notes on every check
  not requiring pandoc, which was unavailable in the environment this was
  verified in - vignette rendering itself was checked separately, by
  running each vignette's extracted R code directly): added the missing
  `fitBatchMix()`-own-formals documentation for
  `auto_tune`/`n_burn`/`mu_proposal_window`/`cov_proposal_window`/
  `m_proposal_window`/`S_proposal_window`/`t_df_proposal_window`
  (`@inheritParams` does not resolve a comma-grouped `@param` tag from
  another function); declared the ggplot2/tidyr non-standard-evaluation
  column names (`.data`, `Acceptance_rate`, `Chain`, `Iteration`,
  `Parameter`, `iteration`, `value`) via `utils::globalVariables()`, the
  standard fix for "no visible binding for global variable" NOTEs on
  NSE-heavy plotting code; added `.Rbuildignore` entries for `.claude`,
  `.Rhistory`, `figure`, `vignettes/figure` and `Rplots.pdf`, and removed
  the (untracked, regenerable) stray copies of the latter two that had
  accumulated locally.
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
