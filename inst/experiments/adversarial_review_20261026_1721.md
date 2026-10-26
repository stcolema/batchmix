# Adversarial review of the samplers and priors, batchmix 3.0.0

Date 2026-10-26. Branch `ccr-5265c537-mfnqub`, from `476bfaa`. Reviewer: Claude (a tool), at the request of Stephen Coleman. Confidence levels are the reviewer's.

## Method

1. Read every C++ sampler (`sampler`, `mvnSampler`, `mvtSampler`, `mvnSamplerSeparationStrategy`, `mvnSamplerMixed`, the four `predictNewBatch*`, `priorOnlyChains`, `genericFunctions`, `pdfs`) and the R code that touches priors, diagnostics and simulation. Each kernel's target, Hastings term and Jacobian was derived by hand and compared with the code.
2. Ran the package's own tests (167 tests, 0 failures before and after the fixes below).
3. Checked the maths by simulation, independently of the package's tests:
   * **Geweke (2004) successive-conditional tests** of the whole sweep (data likelihood included) for MVN, MVT, MVN_LKJ and MVN_MIXED, with sampled `lambda_2` and `rho`, the interaction term and imputed entries, and of the partial-pooling and three Gaussian-process weight updates (`inst/experiments/geweke/`).
   * Finite-difference Jacobian of the correlation reparameterisation, P = 2 to 7.
   * LKJ rejection sampler against the closed-form marginal r ~ 2 Beta(eta - 1 + P/2, eta - 1 + P/2) - 1.
   * `rTruncNorm` against the exact truncated-normal CDF, including the tail switchover at 5 sd.
   * `rankNormalizedRhat` against the `posterior` package.
4. **Not done: the Lean/mathlib confirmation.** The sandbox refused the download and execution of the `elan` installer, and I did not work around that. The algebraic identities I would have formalised (softmax shift invariance; the double-centering projection and its (K-1)(B-1) dimension; conditioning a Gaussian on a sum-to-zero constraint; the exact Gibbs shift mean) were confirmed by hand derivation and by the simulations above instead. That is weaker than a machine-checked proof.

## Defects found (all fixed on this branch, with regression tests in `tests/testthat/test-adversarial-review.R` and `test-predict-new-batch.R`)

| # | Severity | Where | Defect | Evidence | Confidence |
|---|---|---|---|---|---|
| 1 | High | `mvnSamplerMixed` (MVN_MIXED) | `sLogKernel()` summed the batch-scale prior S - S_loc ~ InvGamma over every column. A binary column has S pinned at S_loc, so the prior is evaluated at 0: log(0) and theta/0 give NaN. The acceptance probability of the scale update is then decided by a NaN comparison, and the batch scale of continuous columns was effectively not sampled. | Geweke: successive-conditional mean of S 1.003 against 1.501, acceptance 0.1%. Real fit, three batches with true variance inflation 1, 2.5, 4: acceptance 1-17% and S means 1.01, 1.08, 1.00; after the fix 81-90% and 1.26, 3.84, 4.52 (continuous-only MVN_LKJ reference: 1.26, 4.02, 4.83). After the fix the Geweke test passes (1.502 against 1.501). | High |
| 2 | Medium | `predictNewBatch*.cpp` | The composition sampler is built on the training data pooled with the new batch, so its `delta_2` (empirical-Bayes scale of the shift prior) was that of the pooled data, whereas the fit, and the R code that draws the initial shift, use the training data alone. | By reading; regression test shows the pooled variance is more than twice the fitted one for a shifted batch and that the fitted value is now used. | High |
| 3 | Medium | `generateBatchData(type = "MVT")`, `generateBatchDataMVT()` | One chi-square draw per feature, not per item: independent univariate t margins rather than a multivariate t. | cor(x1^2, x2^2) = 0.001 against 0.15 for a bivariate t with 8 df (0.143 exact). | High |
| 4 | Low | `rStickBreakingPrior()` | Weights sum to 1 - (alpha/(1+alpha))^K on average (0.875 for alpha = 1, K = 3), since the last break was not set to 1 (Ishwaran & James, 2001). Only used to seed initial labels. | Simulation matches the formula. | High |
| 5 | Low | `calcAllocProb()` | Unknown `method` returned `NULL` silently. | By reading. | High |
| 6 | Low | `rLogNormal()` (C++, unused) | Returned `log(randn)`, i.e. NaN for half the draws. Removed. | By reading. | High |

Changing defect 1 changes every fixed-seed MVN_MIXED fit that has a binary column, so `tests/testthat/fixtures/mvn-mixed-semisupervised.rds` was regenerated (only that fixture changed; same names and dimensions, the scale now moves with acceptance 0.55-0.65 and the binary column's scale stays 1). Fits already run with a binary column should be repeated.

## Checked and found correct

* **Metropolis-Hastings corrections**: Gamma random walks for the batch scale, the cluster sds and the t degrees of freedom (forward and reverse densities in the right slots); Wishart random walk for the covariance (including the `-n log|V|` term); symmetric log-scale walks with the +log(value) Jacobian for `rho`, `tau2` and the GP length scale.
* **Priors in the kernels**: the Normal-inverse-Wishart terms (nu + P + 2 exponent), LogNormal sd prior with its 1/sigma factor, LKJ (eta - 1) log det R, the shifted-Gamma prior on the t degrees of freedom, the sum-to-zero Gaussian prior on the interaction and its (K-1)(B-1) degrees of freedom.
* **Reparameterisation of R**: log-Jacobian exact to 1e-9 for P = 2 to 7, so the "empirically fitted" formula is in fact an identity.
* **Empty clusters**: forced prior draws are valid because the order is covariance, then mean (and R, sigma, then mean for LKJ); the joint draw is from the joint prior.
* **Weight priors**: exchangeable softmax logits; the exact Gibbs shift step (mean `-(1/K) sum eta`, covariance Sigma/K for the GP), the sum-to-zero conditioning for `pp_mu` and `gp_beta`, the GLS update of `gp_beta`, and the whitened joint move for the GP hyperparameters (the z-space prior cancels).
* **Imputation**: conditional Gaussian from the precision matrix; the t scale-mixture step (Gamma((nu + P)/2, (nu + d^2)/2) for the mixing variable, then conditional variance 1/(u Lambda_pp)); probit and censored truncation; the Robert (1995) tail sampler.
* **Diagnostics**: Rhat, bulk and tail ESS equal `posterior` to 3e-15 on iid, AR(1), shifted-chain and heavy-tailed inputs.
* The NEWS claim that `pp_mu_prior_sd = 10` corresponds to a Dirichlet concentration of about 0.1: the pairwise log-ratio variance of a Dirichlet(a) is 2 trigamma(a), which equals 200 near a = 0.1.

## Observations (design or documentation; not changed)

1. **Median "allocation probabilities" are not probabilities for K >= 3.** `calcAllocProb()` and `processMCMCChain()` default to the median over draws; the rows of `allocation_probability` then do not sum to 1 (range 0.75-1.0, mean 0.90 in a K = 3 example). `predictNewBatch()` hard-wires the median. The posterior mean is the posterior predictive probability. Changing the default changes results, so it is left to the author.
2. **Data-dependent priors.** `xi`, the inverse-Wishart scale and `delta_2` use the data (empirical Bayes). Reasonable, but it means prior-predictive checks and marginal-likelihood comparisons use the data twice; `kappa = 0.01` (prior sd of a cluster mean ten times the cluster sd) is diffuse enough to make BIC-style comparisons sensitive to it (Lindley-Bartlett effect).
3. The inverse-Wishart default (nu = P + 2) gives the marginal variance an inverse-gamma with shape 1.5, i.e. infinite variance; mean-based checks of it are unreliable (the package's own tests correctly use quantiles).
4. `sampler.h` and `mvnSampler.h` say a uniform rescaling of S is "exactly compensated" by rescaling the cluster covariance. The model multiplies only the diagonal of Sigma_k by S_b, so this holds only for diagonal Sigma_k; fixing the prior mean of S is still a sensible choice, but the stated reason is too strong.
5. `mvnSamplerSeparationStrategy.h` says sigma ~ LogNormal(beta, sqrt(xi)); the code (and kernel) treat `xi` as the sd. They agree at the default xi = 1.
6. `complete_likelihood` is sum of log f(x | theta_label) without the label prior, so it is not the complete-data likelihood; `proposeNewNonNegativeValue()` is unused and has no Hastings correction.

## Limits of this review

* Geweke tests use K = 2, B = 2, P = 2 (K = 3, B = 5 for the weights) and N = 10; larger K or P, the semi-supervised `fixed` path, `continueChain()` and the plotting code were not exercised. The 2020 unit-testing paper cited below reports that Geweke tests can have high false-rejection rates, and consecutive quantile indicators are strongly dependent, so a few moderate z-scores (four beyond 3 among about 800) were not treated as evidence.
* A Geweke test detects disagreement between the sampler and the model as implemented, not an implemented model that differs from the intended one; the latter rests on the hand derivations above.
* One apparent drift (the shift of the first feature in the first batch, positive in 13 of 17 runs) disappeared in four longer runs (10^7 sweeps; z = -2.6, -0.9, 0.2, 0.6) and is treated as chance.

## References checked

* Geweke, J. (2004). Getting it right: joint distribution tests of posterior simulators. *JASA* 99(467), 799-804.
* Ishwaran, H. and James, L. F. (2001). Gibbs sampling methods for stick-breaking priors. *JASA* 96(453), 161-173.
* Lewandowski, D., Kurowicka, D. and Joe, H. (2009). Generating random correlation matrices based on vines and extended onion method. *J. Multivariate Anal.* 100(9), 1989-2001.
* Robert, C. P. (1995). Simulation of truncated normal variables. *Statistics and Computing* 5, 121-125.
* Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C. (2021). Rank-normalization, folding, and localization: an improved R-hat. *Bayesian Analysis* 16(2), 667-718.
* Gelman, A. et al. (2013). *Bayesian Data Analysis*, 3rd ed. (conjugate updates, composition sampling). Stan Reference Manual, "Cholesky factors of correlation matrices" (the transform and its Jacobian). Not re-fetched in this session.
* "Unit testing for MCMC and other Monte Carlo methods", arXiv:2001.06465 (authors not checked), is the source of the false-rejection caveat; I saw only its search-result summary, not the paper, so treat that attribution as unverified.
