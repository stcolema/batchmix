# Regression tests for batchmixControl()/`control = ` (see
# R/batchmixControl.R): bundling the sampler-tuning arguments (proposal
# windows, auto_tune, n_burn) into one object, with the individual
# arguments kept as deprecated, soft-working aliases. Also covers the
# fitBatchMix() `...`-forwarding bug this work surfaced: fitBatchMix()'s
# `...` is documented to forward arbitrary extra arguments on to
# runBatchMix() (include_interaction, r_proposal_window, etc.), but used to
# reject every one of them via an overly strict internal check shared with
# runBatchMix()'s own (deliberately much narrower) `...` contract.

test_that("batchmixControl() validates its arguments and prints", {
  expect_s3_class(batchmixControl(), "batchmix_control")
  expect_output(print(batchmixControl()), "batchmix control")

  expect_error(batchmixControl(mu_proposal_window = "a"), "single numeric")
  expect_error(batchmixControl(mu_proposal_window = c(1, 2)), "single numeric")
  expect_error(batchmixControl(auto_tune = "yes"), "single logical")
  expect_error(batchmixControl(n_burn = -1), "non-negative")
})

test_that("control = is equivalent to the deprecated individual arguments, and the individual arguments still warn", {
  set.seed(1)
  N <- 40
  P <- 2
  X <- matrix(c(rnorm(N * P / 2, 0, 1), rnorm(N * P / 2, 3, 1)), ncol = P, byrow = TRUE)
  batch_vec <- sample(seq(1, 2), replace = TRUE, size = N)
  labels <- sample(0:1, N, replace = TRUE)

  fit_new <- runBatchMix(X, n_iter = 40, thin = 10, batch_vec, "MVN",
    initial_labels = labels, fixed = NULL, verbose = FALSE,
    control = batchmixControl(mu_proposal_window = 0.3, auto_tune = FALSE)
  )
  expect_equal(fit_new$mu_proposal_window, 0.3)
  expect_false(fit_new$auto_tune)

  expect_warning(
    fit_old <- runBatchMix(X, n_iter = 40, thin = 10, batch_vec, "MVN",
      initial_labels = labels, fixed = NULL, verbose = FALSE,
      mu_proposal_window = 0.3, auto_tune = FALSE
    ),
    "deprecated"
  )
  expect_equal(fit_old$mu_proposal_window, 0.3)
  expect_false(fit_old$auto_tune)
})

test_that("when both control and a deprecated individual argument are supplied, control wins (with a warning)", {
  set.seed(2)
  N <- 40
  P <- 2
  X <- matrix(c(rnorm(N * P / 2, 0, 1), rnorm(N * P / 2, 3, 1)), ncol = P, byrow = TRUE)
  batch_vec <- sample(seq(1, 2), replace = TRUE, size = N)
  labels <- sample(0:1, N, replace = TRUE)

  expect_warning(
    fit <- runBatchMix(X, n_iter = 40, thin = 10, batch_vec, "MVN",
      initial_labels = labels, fixed = NULL, verbose = FALSE,
      control = batchmixControl(mu_proposal_window = 0.7), mu_proposal_window = 0.3
    ),
    "using `control`"
  )
  expect_equal(fit$mu_proposal_window, 0.7)
})

test_that("`control` must come from batchmixControl()", {
  set.seed(3)
  N <- 30
  P <- 2
  X <- matrix(c(rnorm(N * P / 2, 0, 1), rnorm(N * P / 2, 3, 1)), ncol = P, byrow = TRUE)
  batch_vec <- sample(seq(1, 2), replace = TRUE, size = N)
  labels <- sample(0:1, N, replace = TRUE)

  expect_error(
    runBatchMix(X, n_iter = 40, thin = 10, batch_vec, "MVN",
      initial_labels = labels, fixed = NULL, verbose = FALSE,
      control = list(mu_proposal_window = 0.3)
    ),
    "batchmixControl"
  )
})

test_that("fitBatchMix()'s `...` forwards arbitrary extra arguments to runBatchMix() instead of erroring", {
  set.seed(4)
  N <- 40
  P <- 2
  X <- matrix(c(rnorm(N * P / 2, 0, 1), rnorm(N * P / 2, 3, 1)), ncol = P, byrow = TRUE)
  batch_vec <- sample(seq(1, 2), replace = TRUE, size = N)
  labels <- sample(0:1, N, replace = TRUE)

  # include_interaction is a formal of runBatchMix(), reached through
  # fitBatchMix()'s `...` - this used to error with "unused argument".
  chains <- fitBatchMix(X, n_chains = 1, n_iter = 40, thin = 10, batch_vec, "MVN",
    initial_labels = labels, fixed = NULL, verbose = FALSE,
    include_interaction = TRUE
  )
  expect_true(chains[[1]]$include_interaction)

  # A genuinely unrecognized argument must still be rejected, just not by
  # fitBatchMix()'s own (too-strict) check - it now surfaces from wherever
  # it is actually detected once forwarded.
  expect_error(
    fitBatchMix(X, n_chains = 1, n_iter = 40, thin = 10, batch_vec, "MVN",
      initial_labels = labels, fixed = NULL, verbose = FALSE,
      totally_bogus_argument = 1
    ),
    "unused argument"
  )
})

test_that("fitBatchMix(): a deprecated MVN_LKJ-only proposal window forwarded via `...` merges into control instead of being silently dropped", {
  set.seed(5)
  N <- 40
  P <- 2
  X <- matrix(c(rnorm(N * P / 2, 0, 1), rnorm(N * P / 2, 3, 1)), ncol = P, byrow = TRUE)
  batch_vec <- sample(seq(1, 2), replace = TRUE, size = N)
  labels <- sample(0:1, N, replace = TRUE)

  expect_warning(
    chains <- fitBatchMix(X, n_chains = 1, n_iter = 40, thin = 10, batch_vec, "MVN_LKJ",
      initial_labels = labels, fixed = NULL, verbose = FALSE,
      r_proposal_window = 0.05
    ),
    "deprecated"
  )
  expect_equal(chains[[1]]$r_proposal_window, 0.05)

  chains_new <- fitBatchMix(X, n_chains = 1, n_iter = 40, thin = 10, batch_vec, "MVN_LKJ",
    initial_labels = labels, fixed = NULL, verbose = FALSE,
    control = batchmixControl(r_proposal_window = 0.05)
  )
  expect_equal(chains_new[[1]]$r_proposal_window, 0.05)
})

test_that("continueChain() resumes proposal windows via control without triggering deprecation warnings", {
  set.seed(6)
  N <- 40
  P <- 2
  X <- matrix(c(rnorm(N * P / 2, 0, 1), rnorm(N * P / 2, 3, 1)), ncol = P, byrow = TRUE)
  batch_vec <- sample(seq(1, 2), replace = TRUE, size = N)
  labels <- sample(0:1, N, replace = TRUE)

  fit1 <- runBatchMix(X, n_iter = 100, thin = 10, batch_vec, "MVN",
    initial_labels = labels, fixed = NULL, verbose = FALSE
  )

  expect_silent(
    fit2 <- continueChain(fit1, X, rep(0L, N), batch_vec, n_iter = 50)
  )
  expect_equal(fit2$n_iter, fit1$n_iter + 50)
})
