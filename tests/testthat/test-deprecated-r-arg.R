test_that("deprecated `R` argument still works and warns, on the main entry points", {
  set.seed(4242)
  N <- 40; P <- 2
  X <- matrix(rnorm(N * P), N, P)
  batch_vec <- sample(1:2, N, replace = TRUE)
  labels <- sample(0:1, N, replace = TRUE)
  fixed <- rep(0, N)

  # The deprecated `R` argument is only supported when either (a) passed
  # positionally in its original slot, or (b) every other argument in the
  # same call is also passed by name - naming `R` while leaving arguments
  # that originally followed it positional is not supported (and cannot be,
  # since `R` no longer occupies that position at all); see `?runBatchMix`.
  set.seed(4242)
  expect_warning(
    out_old <- runBatchMix(X,
      thin = 10, batch_vec = batch_vec, type = "MVN",
      initial_labels = labels, fixed = fixed, R = 50
    ),
    "deprecated"
  )
  set.seed(4242)
  out_new <- runBatchMix(X, 50, thin = 10, batch_vec, "MVN",
    initial_labels = labels, fixed = fixed
  )
  expect_identical(out_old$samples, out_new$samples)
  expect_identical(out_old$n_iter, 50)

  expect_warning(
    fit_old <- batchSemiSupervisedMixtureModel(
      X,
      thin = 10, initial_labels = labels, fixed = fixed, batch_vec = batch_vec,
      type = "MVN", R = 50
    ),
    "deprecated"
  )
  expect_identical(fit_old$n_iter, 50)

  expect_warning(
    chains_old <- fitBatchMix(X, n_chains = 2,
      thin = 10, batch_vec = batch_vec, type = "MVN",
      initial_labels = labels, fixed = fixed, verbose = FALSE, R = 50
    ),
    "deprecated"
  )
  expect_identical(chains_old[[1]]$n_iter, 50)

  expect_warning(
    continued_old <- continueChain(fit_old, X, fixed, batch_vec, R = 20),
    "deprecated"
  )
  expect_identical(continued_old$n_iter, 70)

  expect_warning(
    continued_chains_old <- continueChains(chains_old, X, fixed, batch_vec, R = 20),
    "deprecated"
  )
  expect_identical(continued_chains_old[[1]]$n_iter, 70)

  # Both supplied: n_iter wins, and a (different) warning is still raised.
  expect_warning(
    out_both <- runBatchMix(X, 50, thin = 10, batch_vec, "MVN",
      initial_labels = labels, fixed = fixed, R = 999
    ),
    "n_iter"
  )
  expect_identical(out_both$n_iter, 50)
})

test_that("runMCMCChains() is a deprecated alias for fitBatchMix() and warns", {
  set.seed(99)
  N <- 30; P <- 2
  X <- matrix(rnorm(N * P), N, P)
  batch_vec <- sample(1:2, N, replace = TRUE)
  labels <- sample(0:1, N, replace = TRUE)
  fixed <- rep(0, N)

  set.seed(99)
  expect_warning(
    old_result <- runMCMCChains(X, n_chains = 2, 50, thin = 10, batch_vec, "MVN",
      initial_labels = labels, fixed = fixed, verbose = FALSE
    ),
    "fitBatchMix"
  )
  set.seed(99)
  new_result <- fitBatchMix(X, n_chains = 2, 50, thin = 10, batch_vec, "MVN",
    initial_labels = labels, fixed = fixed, verbose = FALSE
  )
  expect_identical(old_result[[1]]$samples, new_result[[1]]$samples)
})
