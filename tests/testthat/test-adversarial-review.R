# Regression tests for defects found in an adversarial review of the samplers and
# priors. Each was confirmed against an independent check (a Geweke
# successive-conditional run for the sampler, a simulation for the generators)
# before it was fixed.

test_that("MVN_MIXED: the batch scale of a continuous column is updated when another column is binary", {
  # A binary column has its batch scale pinned at S_loc, where the Inverse-Gamma
  # density is not defined (log(0) and theta / 0). sLogKernel() summed that prior
  # over every column, so the kernel was NaN and the scale of every continuous
  # column barely moved (acceptance of a few per cent, posterior means pinned
  # near S_loc however large the batch variance).
  set.seed(5)
  N <- 300
  B <- 3
  batch <- rep(seq_len(B), length.out = N)
  S_true <- c(1, 2.5, 4)
  x1 <- rnorm(N, sd = sqrt(S_true[batch])) + c(0, 0.5, -0.5)[batch]
  X <- cbind(x1, as.numeric(rnorm(N) > 0))
  fit <- runBatchMix(X,
    n_iter = 1500, thin = 1, batch_vec = batch, type = "MVN_MIXED", K_max = 1,
    column_type = c(0, 1), batch_weight_prior = "global", sample_s_scale = FALSE,
    control = batchmixControl(auto_tune = FALSE), verbose = FALSE
  )
  expect_true(all(fit$S_acceptance_rate > 0.4))
  S <- fit$batch_scale[1, , -(1:500)]
  # the noisiest batch has a clearly larger scale than the cleanest one
  expect_gt(mean(S[3, ]), 2 * mean(S[1, ]))
  expect_true(all(fit$batch_scale[2, , ] == 1)) # the binary column's scale stays fixed
})

test_that("rStickBreakingPrior() returns weights that sum to one", {
  set.seed(11)
  for (K in c(1, 3, 10)) {
    w <- replicate(200, rStickBreakingPrior(1, K))
    expect_equal(colSums(matrix(w, nrow = K)), rep(1, 200))
    expect_true(all(w >= 0))
  }
  # the first weight of a Dirichlet-process prior is Beta(1, alpha)
  w1 <- replicate(20000, rStickBreakingPrior(2, 5)[1])
  expect_gt(stats::ks.test(w1, "pbeta", 1, 2)$p.value, 0.001)
})

test_that("generateBatchData(type = 'MVT') shares one chi-square mixing variable across an item's features", {
  # A multivariate t has dependent features even when its scale matrix is diagonal.
  # With W = nu / chi-square, x_j = z_j sqrt(W): E x^2 = nu / (nu - 2), E x1^2 x2^2 =
  # E W^2 = nu^2 / ((nu - 2)(nu - 4)) and E x^4 = 3 E W^2, so for 8 degrees of freedom
  # cor(x1^2, x2^2) = 0.143. Independent t margins give 0.
  set.seed(1)
  dfs <- 8
  nu <- dfs
  expected <- (nu^2 / ((nu - 2) * (nu - 4)) - (nu / (nu - 2))^2) /
    (3 * nu^2 / ((nu - 2) * (nu - 4)) - (nu / (nu - 2))^2)
  d <- generateBatchData(
    N = 60000, P = 2, group_means = 0, group_std_devs = 1, batch_shift = 0, batch_scale = 1,
    group_weights = 1, batch_weights = 1, type = "MVT", group_dfs = dfs, permute_variables = FALSE
  )
  r <- stats::cor(d$observed_data[, 1]^2, d$observed_data[, 2]^2)
  expect_equal(expected, 0.143, tolerance = 0.01)
  expect_gt(r, expected - 0.05)
  expect_lt(r, expected + 0.05)

  d <- generateBatchDataMVT(
    N = 60000, P = 2, group_means = 0, group_std_devs = 1, batch_shift = 0, batch_scale = 1,
    group_weights = 1, batch_weights = 1, dfs = dfs
  )
  r <- stats::cor(d$observed_data[, 1]^2, d$observed_data[, 2]^2)
  expect_gt(r, expected - 0.05)
  expect_lt(r, expected + 0.05)
})

test_that("calcAllocProb() rejects an unknown method instead of returning NULL", {
  chain <- list(
    n_iter = 10, thin = 1, K_max = 2,
    samples = matrix(0L, 10, 3), alloc = array(0.5, c(3, 2, 10))
  )
  expect_error(calcAllocProb(chain, method = "mode"), "must be one of")
})
