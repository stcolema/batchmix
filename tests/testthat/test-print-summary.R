# Regression tests for the S3 print()/summary() methods layered on top of
# the fitting entry points' return values (see R/batchmixFitMethods.R):
# "batchmix_fit" (one chain), "batchmix_fit_list" (several chains) and
# "batchmix_convergence" (assessConvergence()'s return value). These are a
# pure cosmetic addition - the objects are still plain lists underneath, so
# every existing `$`/`[[` access must keep working unchanged; that is
# checked here alongside the new console-facing behaviour.

test_that("runBatchMix()/batchSemiSupervisedMixtureModel() output is classed batchmix_fit, and $ access still works", {
  set.seed(10)
  N <- 40
  P <- 2
  X <- matrix(c(rnorm(N * P / 2, 0, 1), rnorm(N * P / 2, 3, 1)), ncol = P, byrow = TRUE)
  batch_vec <- sample(seq(1, 2), replace = TRUE, size = N)
  labels <- sample(0:1, N, replace = TRUE)

  fit <- runBatchMix(X, n_iter = 60, thin = 10, batch_vec, "MVN",
    initial_labels = labels, fixed = NULL, verbose = FALSE
  )

  expect_s3_class(fit, "batchmix_fit")
  expect_true(is.list(fit))
  expect_identical(fit$type, "MVN")
  expect_equal(fit$N, N)

  expect_output(print(fit), "batchmix fit")
  expect_output(print(fit), "Mean acceptance rates")

  smry <- summary(fit)
  expect_s3_class(smry, "summary.batchmix_fit")
  expect_output(print(smry), "batchmix fit summary")
  expect_output(print(smry), "BIC")
})

test_that("fitBatchMix() output is classed batchmix_fit_list, carries convergence, and prints/summarises", {
  set.seed(11)
  N <- 60
  P <- 2
  X <- matrix(c(rnorm(N * P / 2, 0, 1), rnorm(N * P / 2, 3, 1)), ncol = P, byrow = TRUE)
  batch_vec <- sample(seq(1, 2), replace = TRUE, size = N)
  labels <- sample(0:1, N, replace = TRUE)

  chains <- fitBatchMix(X, n_chains = 3, n_iter = 100, thin = 10, batch_vec, "MVN",
    initial_labels = labels, fixed = NULL, verbose = FALSE
  )

  expect_s3_class(chains, "batchmix_fit_list")
  expect_s3_class(chains[[1]], "batchmix_fit")
  expect_identical(chains[[1]]$type, "MVN")

  expect_output(print(chains), "batchmix fit list")

  smry <- summary(chains)
  expect_s3_class(smry, "summary.batchmix_fit_list")
  expect_length(smry$per_chain, 3)
  expect_output(print(smry), "batchmix fit list summary")

  convergence <- attr(chains, "convergence")
  expect_s3_class(convergence, "batchmix_convergence")
  expect_output(print(convergence), "Convergence")
  expect_match(format(convergence, n_chains = 3), "^3 chains run\\.")
})

test_that("processMCMCChains()/continueChains() keep the batchmix_fit_list/batchmix_fit classes and attributes", {
  set.seed(12)
  N <- 50
  P <- 2
  X <- matrix(c(rnorm(N * P / 2, 0, 1), rnorm(N * P / 2, 3, 1)), ncol = P, byrow = TRUE)
  batch_vec <- sample(seq(1, 2), replace = TRUE, size = N)
  labels <- sample(0:1, N, replace = TRUE)

  chains <- fitBatchMix(X, n_chains = 2, n_iter = 100, thin = 10, batch_vec, "MVN",
    initial_labels = labels, fixed = NULL, verbose = FALSE
  )

  processed <- processMCMCChains(chains, burn = 50)
  expect_s3_class(processed, "batchmix_fit_list")
  expect_s3_class(processed[[1]], "batchmix_fit")
  expect_identical(attr(processed, "best_chain"), attr(chains, "best_chain"))

  extended <- continueChains(chains, X, rep(0L, N), batch_vec, n_iter = 50)
  expect_s3_class(extended, "batchmix_fit_list")
  expect_s3_class(extended[[1]], "batchmix_fit")
})

test_that("getBestChain() returns a single batchmix_fit with its convergence attached", {
  set.seed(13)
  N <- 60
  P <- 2
  X <- matrix(c(rnorm(N * P / 2, 0, 1), rnorm(N * P / 2, 3, 1)), ncol = P, byrow = TRUE)
  batch_vec <- sample(seq(1, 2), replace = TRUE, size = N)
  labels <- sample(0:1, N, replace = TRUE)

  chains <- fitBatchMix(X, n_chains = 3, n_iter = 100, thin = 10, batch_vec, "MVN",
    initial_labels = labels, fixed = NULL, verbose = FALSE
  )

  best <- getBestChain(chains)
  expect_s3_class(best, "batchmix_fit")
  expect_s3_class(attr(best, "convergence"), "batchmix_convergence")
  expect_output(print(best), "Convergence")
})
