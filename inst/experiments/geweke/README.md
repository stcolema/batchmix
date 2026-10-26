# Geweke successive-conditional tests of the samplers

Geweke (2004, "Getting it right: joint distribution tests of posterior simulators",
JASA 99(467)) compares two ways of drawing from the joint law of parameters and data.

* **Marginal-conditional**: draw the parameters from the prior, then the data given the
  parameters. Independent draws.
* **Successive-conditional**: alternate "draw data given parameters" with one sweep of the
  sampler (which targets the parameters given data). If every kernel is correct, the
  chain's stationary law of the parameters is the prior, so the two must agree.

Unlike the prior-only chains in `src/priorOnlyChains.cpp`, the data likelihood stays in
every kernel, so this tests the likelihood, the Hastings corrections, the cached
matrices, the empty-cluster draws and the label and weight updates together. The
data-dependent hyperparameters (`xi`, `scale`, `delta_2`, ...) are fixed by hand, since
they would otherwise move with the simulated data.

| Script (driver, source) | What it runs |
|---|---|
| `run_geweke.R`, `geweke_mvn.cpp` | MVN and MVT, K = 2, B = 2, P = 2, N = 10; optionally sampled `lambda_2` and `rho`, the batch-by-cluster interaction, and imputed (missing) entries |
| `run_geweke_sep.R`, `geweke_sep.cpp` | MVN_LKJ, and MVN_MIXED (one continuous column left-censored at -0.7, one probit column) |
| `run_geweke_w.R`, `geweke_w.cpp` | the weight update alone: partial pooling and the Matern-3/2, RW1 and RW2 Gaussian processes, with sampled hyperparameters, K = 3, B = 5 |

Run from this directory (the scripts compile with `Rcpp::sourceCpp()` against `../../../src`):

```sh
Rscript run_geweke.R <is_t> <n_sc> <thin> <sample_hyper> <sample_hyper> <seed> <interaction> <missing>
Rscript run_geweke_sep.R <mixed> <n_sc> <thin> <seed>
Rscript run_geweke_w.R <type: 1 pp, 2 gp> <kernel: 0, 1, 2> <sample_hyper> <n_sc> <seed>
```

For each summary the output gives the Monte Carlo and successive-conditional means with a
z-score (batch-means standard error for the successive-conditional chain), and z-scores of
the indicators `1{x <= q}` at the marginal-conditional 5, 25, 50, 75 and 95% quantiles.
Read the quantile columns first: several priors here (the inverse-Wishart diagonal) have
infinite variance, so mean comparisons for `cov0_00` are unreliable. The `N0` row (number
of items in cluster 0) is discrete, so its quantile z-scores are meaningless.

## Result when this was written (2026-10-26)

Fifteen runs (MVN, MVT, MVN_LKJ, MVN_MIXED, and MVN/MVT with the interaction or with
imputed entries; 4 x 10^5 marginal-conditional draws and 1.5 x 10^6 to 6 x 10^6 sweeps
each) and five weight-update runs.

* **One discrepancy, now fixed.** MVN_MIXED with a binary column: the batch scale of the
  continuous column sat near its lower bound in the successive-conditional chain (mean 1.003
  against 1.501 under the prior; acceptance rate of the update 0.1%), because the kernel
  was NaN (the Inverse-Gamma prior evaluated at the binary column's pinned scale). After
  the fix: 1.502 against 1.501, acceptance 84%.
* **Everything else agrees.** About 800 z-scores outside the `N0` rows: largest |z| of a
  mean 2.6, of a quantile indicator 3.5, four of them above 3 where independence would give
  about two. The indicators at neighbouring quantiles are strongly positively dependent
  and the batch-means standard error from 40 batches is itself noisy, so a small excess of
  moderate z-scores is expected and is not evidence of a defect.
* **One pattern, checked and not confirmed.** The batch shift of the first feature in the
  first batch (`m00`) had a positive mean z in 13 of 17 earlier runs (pooled about 2.8
  standard errors, chosen after looking at about a dozen statistics). Four further MVN runs
  of 10^7 sweeps each gave z = -2.6, -0.9, 0.2 and 0.6 (mean -0.2): the pattern went away.
  The shift of the second feature in the second batch never showed it. Treat it as chance.
* The weight-update runs (partial pooling, Matern-3/2, RW1, RW2) show no discrepancy
  (largest |z| of a mean 1.8, of a quantile indicator 3.0).
