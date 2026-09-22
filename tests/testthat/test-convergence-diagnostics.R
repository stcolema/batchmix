# Regression tests for rankNormalizedRhat()/assessConvergence(), pinning
# the specific behaviour that motivates using the rank-normalized, folded
# diagnostic (Vehtari et al., 2021) over the classical Gelman-Rubin Rhat:
# it must flag both location AND scale disagreement across chains, even
# when the other of the two agrees.

test_that("rankNormalizedRhat: well-mixed chains report Rhat close to 1 and full ESS", {
  set.seed(1)
  good_chains <- matrix(rnorm(3000), ncol = 3)
  out <- rankNormalizedRhat(good_chains)

  expect_lt(abs(out$rhat - 1), 0.02)
  expect_gt(out$ess_bulk, 2000)
  expect_gt(out$ess_tail, 2000)
})

test_that("rankNormalizedRhat: chains with different means are flagged (bulk)", {
  set.seed(2)
  bad_chains <- cbind(rnorm(1000, 0), rnorm(1000, 2), rnorm(1000, -2))
  out <- rankNormalizedRhat(bad_chains)

  expect_gt(out$rhat, 1.1)
  expect_gt(out$rhat_bulk, 1.1)
})

test_that("rankNormalizedRhat: chains with the same mean but different scales are flagged (tail), not by classical bulk Rhat", {
  set.seed(3)
  scale_chains <- cbind(rnorm(1000, 0, 1), rnorm(1000, 0, 5), rnorm(1000, 0, 0.2))
  out <- rankNormalizedRhat(scale_chains)

  # The whole point of folding: bulk agrees (same location) but tail must not.
  expect_lt(out$rhat_bulk, 1.05)
  expect_gt(out$rhat_tail, 1.1)
  expect_equal(out$rhat, max(out$rhat_bulk, out$rhat_tail))
})

test_that("rankNormalizedRhat accepts a list of unequal-length chains", {
  set.seed(4)
  chains <- list(rnorm(500), rnorm(520), rnorm(480))
  out <- rankNormalizedRhat(chains)

  expect_true(is.finite(out$rhat))
  expect_gt(out$rhat, 0)
})

test_that("assessConvergence + getBestChain integrate with fitBatchMix", {
  set.seed(5)
  N <- 60
  P <- 2
  X <- matrix(c(rnorm(N * P / 2, 0, 1), rnorm(N * P / 2, 3, 1)), ncol = P, byrow = TRUE)
  batch_vec <- sample(seq(1, 2), replace = TRUE, size = N)
  labels <- sample(0:1, N, replace = TRUE)

  chains <- fitBatchMix(X, n_chains = 3, n_iter = 200, thin = 20, batch_vec, "MVN",
    initial_labels = labels, fixed = NULL, verbose = FALSE
  )

  convergence <- attr(chains, "convergence")
  expect_false(is.null(convergence))
  expect_true(is.finite(convergence$rhat))
  expect_true(convergence$best_chain %in% seq_len(3))
  expect_identical(attr(chains, "best_chain"), convergence$best_chain)

  best <- getBestChain(chains)
  expect_identical(best$samples, chains[[convergence$best_chain]]$samples)

  # processMCMCChains carries the attributes through; continueChains
  # recomputes them on the extended chains.
  processed <- processMCMCChains(chains, burn = 50)
  expect_identical(attr(processed, "best_chain"), attr(chains, "best_chain"))

  extended <- continueChains(chains, X, rep(0L, N), batch_vec, n_iter = 100)
  expect_false(is.null(attr(extended, "best_chain")))
})
