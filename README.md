# batchmix

Downloads from CRAN in the past month:

[![](https://cranlogs.r-pkg.org/badges/batchmix)](https://cran.r-project.org/package=batchmix)

Semi-supervised and unsupervised Bayesian mixture models that simultaneously
infer the cluster/class structure and a batch correction. Densities available
are the multivariate normal (with an inverse-Wishart or an LKJ prior on the
cluster covariance), the multivariate t, and a mixed continuous/binary/censored
model with missing data. The sampler is implemented in C++. This package is
aimed at analysis of low-dimensional data generated across several batches. See
[Coleman et al. (2022)](https://doi.org/10.1101/2022.01.14.476352) for details
of the model.

Installation needs a C++ compiler (as for any package with compiled code) and
no other system tooling: there is no Rust dependency.

## Read this if you are updating from an earlier version

This release changes defaults and results, so analyses run with an earlier
version will not reproduce. The same notice is printed by `library(batchmix)`.

* **Partial pooling is on by default.** Batches get their own mixture weights,
  drawn exchangeably around an estimated population composition
  (`batch_weight_prior = "partial_pooling"`; `"global"` when there is a single
  batch), and the concentration of the batch-scale prior is estimated
  (`sample_s_scale = TRUE`). For the earlier fixed-hyperparameter model use
  `batch_weight_prior = "global", sample_s_scale = FALSE`. Under partial
  pooling `alpha`/`concentration` are unused, so a generous `K_max` is no
  longer pushed towards empty components by a sparse Dirichlet prior; choose
  `K_max` near the number of clusters you expect.
* **`batch_corrected_data` / `inferred_dataset`** are now the posterior mean of
  each item's batch-free signal, the exact inverse of the fitted batch model.
  They are shrunk towards the item's cluster mean (a denoised estimate), not a
  rescaled residual.
* **Sampler corrections change fixed-seed output.** Proposal-window
  auto-tuning previously moved the covariance, batch-scale, marginal-scale and
  degrees-of-freedom windows in the wrong direction; the `MVN_LKJ` /
  `MVN_MIXED` prior on the marginal standard deviations was mis-specified;
  `mu_proposal_window` is now in units of each cluster's spread. See
  `NEWS.md`.
* **Compare models with `calcBICM()`** (BIC at the best retained draw), not the
  per-iteration `BIC` trace, which uses a single posterior draw.
* **`salso` is no longer a dependency.** Unsupervised point estimates use
  `minVI()` on the posterior similarity matrix, which needs the `N x N` matrix
  in memory (`8 N^2` bytes); `nCores` is ignored.
* `fitBatchMix()` is the main entry point (`runMCMCChains()` is a deprecated
  alias) and the iteration count `R` is now `n_iter`.

## Where to start

The vignettes are the documentation of the workflow; please read them before
re-running an old analysis (`browseVignettes("batchmix")`):

| Vignette | What it covers |
|---|---|
| `batchmix_workflow` | An end-to-end analysis: simulate, prior predictive check, fit several chains, diagnose, posterior predictive check, recover the truth. **Start here.** |
| `batch_weight_priors` | What partial pooling changes: global, partially pooled and Gaussian-process priors on the class proportions, and pooling of the batch scales. |
| `covariance_models` | Choosing between `MVN`, `MVT` and `MVN_LKJ`, comparing numbers of clusters with BICM, label switching. |
| `probit_missing_censored` | Mixed continuous/binary data with missing and censored entries (`MVN_MIXED`). |

A minimal analysis:

```r
library(batchmix)

# X: N x P matrix; batch_vec: batch label per row
fit <- fitBatchMix(X, n_chains = 4, n_iter = 5000, thin = 25,
                   batch_vec = batch_vec, type = "MVN", K_max = 3)

fit                                   # Rhat / ESS, best chain
best <- getBestChain(fit)
calcBICM(best)                        # for comparing models
est <- processMCMCChain(best, burn = 2500)
est$inferred_dataset                  # batch-corrected data
```

## Advice on using the package

* **Run several chains and check convergence.** `fitBatchMix()` reports
  rank-normalized split-Rhat and effective sample size on the complete-data
  log-likelihood (which is invariant to label switching). If Rhat is above about
  1.01, continue the chains with `continueChains()` or run more of them; do
  not read a point estimate off chains that disagree.
* **Proposal windows are tuned automatically.** Metropolis-Hastings proposal
  windows are adapted during burn-in (Robbins-Monro, targeting an acceptance
  rate of 0.234 for block updates and 0.44 for the degrees of freedom); the
  windows and `auto_tune`/`n_burn` live in `batchmixControl()`. Check the
  acceptance rates afterwards with `plotAcceptanceRates()`; rates between about
  0.1 and 0.5 are healthy. Empty clusters are drawn from their prior and are
  not counted as acceptances.
* **Look at the priors before fitting** with `simulatePriorPredictive()` (it
  simulates the default partial-pooling model) and at the fit afterwards with
  `simulatePosteriorPredictive()` / `plotPredictiveCheck()`.
* **Missing data are handled inside the sampler** (mark them `NA`), and
  `type = "MVN_MIXED"` adds probit-linked binary columns and censored entries.
  The likelihood behind `BIC`/`calcBICM()` includes the imputed/latent values, so
  do not compare those criteria between models that treat missing, censored or
  binary entries differently.
* **New batches.** `predictNewBatch()` draws the parameters and class
  probabilities of an unseen batch conditional on the fitted model.
