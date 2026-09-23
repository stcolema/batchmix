#!/usr/bin/Rscript
# Regression test for a crash reported against type = "MVN_LKJ" (and, by
# inheritance, "MVN_MIXED"): `inv_sympd(): matrix is singular or not
# positive definite`, later `Mat::operator(): index out of bounds`.
#
# Forensic summary (see NEWS.md for the full account): a random-walk
# proposal for the correlation matrix R, reparameterised via a Cholesky/
# partial-correlation transform, is guaranteed PD in exact arithmetic for
# any candidate in (-1,1)^(P(P-1)/2) - but that does not keep it away from
# the boundary of the PD cone, and a proposal legitimately close to that
# boundary can be PD in exact arithmetic yet numerically singular in
# floating point (confirmed via gdb backtrace on a real run: eigenvalues
# down to ~1.4e-4, with the smallest nominally -1.3e-16 after
# symmetrising). Two separate, compounding bugs followed from this:
#   1. Every `inv_sympd(proposed_cov)`-style call across
#      mvnSamplerSeparationStrategy.cpp/mvnSamplerMixed.cpp's Metropolis
#      steps (rMHStep, sigmaMHStep, batchScaleMetropolis) and
#      sampleCovPrior()'s initial draw used the throwing single-argument
#      form, uncaught - a numerically-degenerate proposal crashed the
#      whole chain instead of simply being (or needing to be) rejected.
#   2. rMHStep()'s regular (non-empty-cluster) branch called the
#      non-throwing, output-parameter form of arma::chol() but never
#      checked its boolean return value; on failure that form leaves the
#      output EMPTY (0x0) rather than P x P, and the very next line
#      (choleskyToPartialCorrelations()) indexed it assuming P x P -
#      an out-of-bounds crash, not a hypothetical (confirmed via gdb).
#
# This test reproduces the exact scenario that crashed before both fixes
# (four chains, four correlated/uncorrelated features, unsupervised
# MVN_LKJ) - it is necessarily slowish (a genuine, seed-dependent
# floating-point rare event, not a deterministic one-line repro), but a
# smaller/faster version was not found to reproduce the crash reliably.
test_that("fitBatchMix(..., type = 'MVN_LKJ') does not crash on a scenario that previously threw inv_sympd()/index-out-of-bounds errors", {
  skip_on_cran()

  set.seed(303)
  N <- 500
  P <- 4
  K <- 2
  B <- 2

  R_true <- list(
    matrix(c(
      1, 0.7, 0, 0,
      0.7, 1, 0, 0,
      0, 0, 1, -0.5,
      0, 0, -0.5, 1
    ), 4, 4),
    diag(4)
  )
  sigma_true <- rbind(c(1, 1, 1, 1), c(1.2, 1.2, 1.2, 1.2))
  mu_true <- rbind(c(0, 0, 0, 0), c(4, 4, 3, 3))
  batch_shift_true <- rbind(c(0, 0, 0, 0), c(0.4, -0.3, 0.2, 0.1))

  labels <- sample(0:(K - 1), N, replace = TRUE)
  batch_vec <- sample(0:(B - 1), N, replace = TRUE)

  X <- matrix(0, N, P)
  for (i in seq_len(N)) {
    k <- labels[i] + 1
    b <- batch_vec[i] + 1
    Sigma_k <- diag(sigma_true[k, ]) %*% R_true[[k]] %*% diag(sigma_true[k, ])
    X[i, ] <- mu_true[k, ] + batch_shift_true[b, ] +
      as.numeric(MASS::mvrnorm(1, rep(0, P), Sigma_k))
  }

  expect_no_error(
    fit_lkj <- fitBatchMix(
      X, n_chains = 4, n_iter = 4000, thin = 20, batch_vec, "MVN_LKJ",
      K_max = K, initial_labels = labels, fixed = NULL, eta = 1.0
    )
  )

  expect_length(fit_lkj, 4)
  expect_true(all(vapply(fit_lkj, function(ch) all(is.finite(ch$BIC)), logical(1))))
})
