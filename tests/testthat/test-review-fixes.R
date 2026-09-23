#!/usr/bin/Rscript
# Regression tests for a batch of correctness issues found during a
# mathematical/Bayesian-rigour review of the package (2026-09) and fixed in
# the same change: the semi-supervised observed-data likelihood/BIC not
# conditioning on known (fixed) labels, minVI()'s max.k crash, the
# mean_sum_est/cov_comb_est batch-vs-cluster column indexing bug in
# processMCMCChain(), the burn-in off-by-one, prepareInitialParameters()'s
# `&`/`|` dimension-check bug, relabelChain() - the label-switching
# correction that processMCMCChain()'s point estimates previously lacked
# entirely - and, once relabelChain() made it safe to do so,
# processMCMCChain() no longer withholding pred/prob/allocation_probability
# for unsupervised fits.

test_that("fixed (semi-supervised) items' observed_likelihood matches complete_likelihood + log w, not a marginalisation over every component (regression for the fixed-label BIC inflation bug)", {
  # When every item's label is fixed, updateAllocation() never marginalises
  # over components for any item, so observed_likelihood and
  # complete_likelihood become a deterministic function of (labels, ll, w)
  # at every iteration: observed_likelihood(t) - complete_likelihood(t) must
  # equal sum_n log(w(labels_n)) exactly (up to floating-point tolerance).
  # Before the fix, observed_likelihood instead log-sum-exp'd over every
  # component regardless of `fixed`, which is strictly larger and would
  # break this identity.
  set.seed(20260101)
  N <- 80; P <- 1; K <- 2; B <- 2
  X <- matrix(c(rnorm(N / 2, 0, 1), rnorm(N / 2, 8, 1)), ncol = 1)
  labels <- c(rep(0, N / 2), rep(1, N / 2))
  batch_vec <- sample(0:(B - 1), N, replace = TRUE)

  out <- runBatchMix(X, 200, 10, batch_vec, "MVN",
    initial_labels = labels, fixed = rep(1, N)
  )

  n_saved <- length(out$observed_likelihood)
  # weights(t, k+1) is the sampled weight for (0-indexed) component k at
  # iteration t; labels never change for fully-fixed items, so `labels` is
  # simply the initial vector throughout.
  expected_log_w_sum <- vapply(seq_len(n_saved), function(t) {
    sum(log(out$weights[t, labels + 1]))
  }, numeric(1))

  actual_gap <- out$observed_likelihood[, 1] - out$complete_likelihood[, 1]

  expect_equal(actual_gap, expected_log_w_sum, tolerance = 1e-6)
})

test_that("minVI() does not crash when max.k is passed explicitly (regression for the undefined k_inds bug)", {
  set.seed(1)
  n <- 20
  cl <- rep(1:4, each = n / 4)
  psm <- createSimilarityMat(matrix(rep(cl, 10), nrow = 10, byrow = TRUE))

  expect_no_error(res_avg <- minVI(psm, method = "avg", max.k = 6))
  expect_no_error(res_comp <- minVI(psm, method = "comp", max.k = 6))
  expect_length(res_avg, n)
  expect_length(res_comp, n)
  expect_equal(attr(res_avg, "info")$maxNClusters, 6)
})

test_that("processMCMCChain()'s mean_sum_est/cov_comb_est index (cluster, batch) pairs correctly for K_max > 1 and B > 1 (regression for the contiguous-vs-strided column indexing bug)", {
  # mean_sum/cov_comb are stored with columns for (cluster k, batch b) at
  # position (k-1)*B + b - i.e. strided by B for a fixed batch, not
  # contiguous. A wrong (contiguous) extraction would pull data for the
  # wrong (cluster, batch) pair whenever K_max > 1 and B > 1, so
  # mean_sum_est[, k, b] would no longer equal mean_est[, k] + shift_est[, b]
  # (the very definition of mean_sum = mu_k + m_b).
  set.seed(20260102)
  P <- 1; K <- 3; B <- 3
  mu_true <- c(-6, 0, 6); m_true <- c(-2, 0, 2)

  X <- c(); labels <- c(); batch_vec <- c()
  n_per_cell <- 40
  for (k in 1:K) for (b in 1:B) {
    X <- c(X, rnorm(n_per_cell, mu_true[k] + m_true[b], 0.3))
    labels <- c(labels, rep(k - 1, n_per_cell))
    batch_vec <- c(batch_vec, rep(b - 1, n_per_cell))
  }
  X <- matrix(X, ncol = 1)

  out <- runBatchMix(X, 500, 10, batch_vec, "MVN",
    initial_labels = labels, fixed = rep(1, length(labels))
  )

  # mean_sum(t) = mu(t) + m(t) exactly, every iteration, by construction in
  # the C++ sampler - but only the MEAN is linear (mean(mu) + mean(m) ==
  # mean(mu + m) exactly); the (default) median is not, so use "mean" here
  # to make the identity an exact, tight check of the (cluster, batch)
  # indexing rather than a loose one confounded by median's non-linearity.
  proc <- processMCMCChain(out, burn = 200, point_estimate_method = "mean")

  for (k in seq_len(K)) {
    for (b in seq_len(B)) {
      expect_equal(
        proc$mean_sum_est[1, k, b],
        proc$mean_est[1, k] + proc$shift_est[1, b],
        tolerance = 1e-8,
        label = sprintf("mean_sum_est[1, %d, %d]", k, b)
      )
      # cov_comb(p,p) = cov(p,p) * S(p,b) exactly every iteration, but
      # mean(cov * S) != mean(cov) * mean(S) in general (a product, not a
      # sum) - so this is a looser, approximate check than the mean_sum one
      # above (E[XY] ~= E[X]E[Y] here only because Cov(cov, S) is small over
      # a well-mixed chain). Still tight enough to fail hard if the
      # (cluster, batch) indexing were wrong: a scrambled block would pull
      # an entirely different cluster's covariance and/or batch's scale,
      # not a value a few percent off.
      expect_equal(
        proc$cov_comb_est[[b]][1, 1, k],
        proc$cov_est[1, 1, k] * proc$scale_est[1, b],
        tolerance = 0.1,
        label = sprintf("cov_comb_est[[%d]][1,1,%d]", b, k)
      )
    }
  }
})

test_that("processMCMCChain() applies no burn-in when 0 < burn < thin (regression for seq(1, 0) == c(1, 0))", {
  set.seed(20260103)
  N <- 40; P <- 1
  X <- matrix(rnorm(N), ncol = 1)
  labels <- sample(0:1, N, replace = TRUE)
  batch_vec <- rep(0, N)

  out <- runBatchMix(X, 500, 50, batch_vec, "MVN",
    initial_labels = labels, fixed = rep(0, N)
  )

  n_saved_raw <- dim(out$means)[3]
  # burn = 10 < thin = 50, so floor(burn / thin) == 0: nothing should be
  # dropped. seq(1, 0) == c(1, 0) would have dropped the first saved sample.
  proc <- processMCMCChain(out, burn = 10)

  expect_equal(dim(proc$means)[3], n_saved_raw)
  expect_equal(nrow(proc$samples), n_saved_raw)
})

test_that("prepareInitialParameters() rejects an initial matrix with only one dimension wrong (regression for the `&` instead of `|` bug)", {
  P <- 2; K <- 3; B <- 2

  # Right number of columns (K), wrong number of rows (P) - the old `&`
  # check required BOTH dimensions to be wrong to trigger, so this slipped
  # through silently before.
  wrong_rows <- matrix(0, nrow = P + 1, ncol = K)
  expect_error(
    prepareInitialParameters(wrong_rows, NULL, NULL, NULL, NULL, P, K, B, "MVN"),
    "Initial class means"
  )

  wrong_shift_rows <- matrix(0, nrow = P + 1, ncol = B)
  expect_error(
    prepareInitialParameters(NULL, NULL, wrong_shift_rows, NULL, NULL, P, K, B, "MVN"),
    "Initial batch shifts"
  )

  wrong_scale_cols <- matrix(1, nrow = P, ncol = B + 1)
  expect_error(
    prepareInitialParameters(NULL, NULL, NULL, wrong_scale_cols, NULL, P, K, B, "MVN"),
    "Initial batch scales"
  )

  # A correctly-shaped matrix should still be accepted.
  ok <- matrix(0, nrow = P, ncol = K)
  expect_no_error(prepareInitialParameters(ok, NULL, NULL, NULL, NULL, P, K, B, "MVN"))
})

test_that("relabelChain() aligns cluster-indexed arrays across a synthetic label switch (the centrepiece fix for label-switching-unsafe point estimates)", {
  # Two "iterations" of a K=2, P=1, B=1 chain, N=4 items, where the
  # component roles are exactly swapped between iterations - a clean,
  # deterministic stand-in for genuine MCMC label switching. Before this
  # fix, processMCMCChain() averaged such arrays directly by raw component
  # index, which for real chains blurs two genuinely different components
  # towards their midpoint instead of recovering either one.
  #
  # `samples` is 0-INDEXED (0..K_max-1), matching the raw C++ sampler's own
  # convention (mcmc_output$samples in a real fit is never 1-indexed) -
  # this is deliberately not 1, 2 here, since a synthetic 1-indexed fixture
  # was exactly what let a real 0-vs-1-indexed bug (samples containing a 0
  # made `perm[samples[t, ]]` drop that element via R's "index 0" rule,
  # confirmed via a real R CMD check failure on plotBatchCorrection()'s own
  # example) slip past this same test previously.
  samples <- rbind(
    c(0, 0, 1, 1), # iteration 1: items 1-2 -> component 0, items 3-4 -> component 1
    c(1, 1, 0, 0)  # iteration 2 (reference): the same grouping, swapped labels
  )

  means <- array(0, dim = c(1, 2, 2))
  means[, , 1] <- c(10, 20) # iter 1: component 0 = 10, component 1 = 20
  means[, , 2] <- c(20, 10) # iter 2: component 0 = 20, component 1 = 10 (same clusters, swapped labels)

  mcmc_output <- list(samples = samples, means = means, K_max = 2)

  relabelled <- relabelChain(mcmc_output)

  # Aligned to the reference (iteration 2): both iterations should now
  # agree exactly that component 1 (R array position, i.e. label 0) = 20
  # and component 2 (label 1) = 10.
  expect_equal(relabelled$means[, , 1], c(20, 10))
  expect_equal(relabelled$means[, , 2], c(20, 10))
  expect_equal(relabelled$samples[1, ], c(1, 1, 0, 0))
  expect_equal(relabelled$samples[2, ], c(1, 1, 0, 0))

  # The naive (unaligned) average would have blurred both components to
  # their midpoint, 15 - masking the true bimodal structure entirely.
  naive_mean <- rowMeans(means, dims = 2L)
  expect_true(all(abs(naive_mean - 15) < 1e-8))

  # The relabelled average recovers the true, well-separated component
  # values exactly (this toy example is noiseless).
  aligned_mean <- rowMeans(relabelled$means, dims = 2L)
  expect_equal(sort(as.vector(aligned_mean)), c(10, 20))
})

test_that("relabelChain() leaves already-aligned chains and single-component/single-iteration input unchanged", {
  samples <- rbind(c(0, 0, 1, 1), c(0, 0, 1, 1))
  means <- array(0, dim = c(1, 2, 2))
  means[, , 1] <- c(5, -5)
  means[, , 2] <- c(5, -5)

  out <- relabelChain(list(samples = samples, means = means, K_max = 2))
  expect_equal(out$means, means)
  expect_equal(out$samples, samples)

  # K_max < 2: nothing to relabel, input returned as-is.
  one_k <- list(samples = matrix(0, 2, 3), means = array(1, dim = c(1, 1, 2)), K_max = 1)
  expect_identical(relabelChain(one_k), one_k)
})

test_that("relabelChain()/processMCMCChain() work on a real MVN fit's raw (0-indexed) output (regression for two compounding bugs found via R CMD check --run-donttest on plotBatchCorrection()'s own example)", {
  # Two distinct bugs, both only visible on REAL sampler output (every
  # synthetic fixture above used 1-indexed `samples`, which is why neither
  # was caught until a real R CMD check run):
  #  1. mcmc_output$samples is 0-indexed (0..K_max-1) - relabelChain() was
  #     written and tested assuming 1-indexed labels, so `perm[samples[t,]]`
  #     silently dropped any item with label 0 (R indexing with 0 selects
  #     nothing), corrupting `samples`' length.
  #  2. `mcmc_output$t_df` (accessed via `$`, which does partial matching)
  #     silently matched `t_df_proposal_window` for a non-MVT fit instead
  #     of returning NULL, since there is no field literally named "t_df"
  #     on an MVN fit - `has_t_df` was wrongly TRUE, and every downstream
  #     `t_df` access then operated on a bare scalar instead of an
  #     (n_saved x K) matrix.
  set.seed(1)
  X <- matrix(c(rnorm(100, 0, 1), rnorm(100, 3, 1)), ncol = 2, byrow = TRUE)
  batch_vec <- sample(seq(1, 3), replace = TRUE, size = 100)

  mcmc_out <- runBatchMix(X, 1000, 50, batch_vec, "MVN", verbose = FALSE)
  expect_equal(min(mcmc_out$samples), 0) # sanity: genuinely 0-indexed
  expect_lt(max(mcmc_out$samples), mcmc_out$K_max)

  expect_no_error(proc <- processMCMCChain(mcmc_out, burn = 250))
  expect_equal(dim(proc$samples), dim(mcmc_out$samples[-seq_len(5), , drop = FALSE]))

  expect_no_error(p <- plotBatchCorrection(mcmc_out, X, batch_vec, burn = 250))
  expect_s3_class(p, "ggplot")
})

test_that("processMCMCChain() computes allocation_probability/prob/pred for a fully unsupervised fit, not only semi-supervised ones (regression for an unwarranted is_semisupervised gate)", {
  # Before relabelChain() existed, gating pred/prob on Semisupervised may
  # have been a deliberate (if undocumented) safety measure, since a fully
  # unsupervised fit's raw allocation draws are not label-switching-safe to
  # average. relabelChain() removes that reason - alloc is relabelled by
  # processMCMCChain() before calcAllocProb() is called - so this gate was
  # left over, unconditionally withholding pred/prob for any unsupervised
  # fit (e.g. exactly the "does the clustering recover known labels"
  # comparison the covariance-model vignette makes for MVN vs MVT).
  set.seed(20260301)
  N <- 150; P <- 1
  true_means <- c(0, 8)
  labels <- sample(0:1, N, replace = TRUE)
  X <- matrix(rnorm(N, true_means[labels + 1], 0.5), ncol = 1)
  batch_vec <- rep(0, N)

  out <- runBatchMix(X, 400, 10, batch_vec, "MVN", K_max = 2)
  expect_false(out$Semisupervised)

  proc <- processMCMCChain(out, burn = 200)

  expect_false(is.null(proc$pred))
  expect_length(proc$pred, N)
  expect_equal(dim(proc$allocation_probability), c(N, 2))
  expect_length(proc$prob, N)

  # A real (if weak) ground-truth check: predicted labels should align with
  # the true, well-separated generating classes far better than chance,
  # accounting for the arbitrary label permutation.
  tab <- table(proc$pred, labels)
  accuracy <- max(sum(diag(tab)), sum(tab) - sum(diag(tab))) / sum(tab)
  expect_gt(accuracy, 0.85)
})
