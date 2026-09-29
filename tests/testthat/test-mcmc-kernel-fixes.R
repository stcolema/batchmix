# Regression tests for four sampler defects found in a code review:
#   1. Robbins-Monro adaptation moved the reciprocal-parameterised windows
#      (cov, S, t_df, sigma) the wrong way.
#   2. sigmaLogKernel() omitted the log-normal Jacobian for its Gamma random
#      walk on sigma, so the sampler targeted LogNormal(beta + xi^2, xi)
#      rather than LogNormal(beta, xi).
#   3. r_log_det was never initialised, so the first R Metropolis decisions
#      scored the current state with log det(R) = 0.
#   4. Empty-cluster (forced prior) draws were counted as acceptances, and
#      the MVN_LKJ empty-cluster draw of mu preceded the draw of its own
#      covariance, so (mu, Sigma) did not follow the joint prior.
#
# The prior-only chains use priorOnlyLKJChain(), which swaps the data
# likelihood for zero so each chain's stationary law is the prior; the
# expected marginals are known in closed form.

lkj_prior_chain <- function(K = 2, B = 2, n_iter = 20000, eta = 1.0, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(10), 5, 2)
  priorOnlyLKJChain(
    X = X, K = K, B = B,
    labels = rep(0L, 5), # cluster 0 occupied, every other cluster empty
    batch_vec = c(0L, 1L, 0L, 1L, 0L)[seq_len(5)] %% B,
    n_iter = n_iter, eta = eta,
    r_pw = 1.0, sigma_pw = 5, mu_pw = 8, m_pw = 0.3, S_pw = 5
  )
}

test_that("sigma Metropolis step targets LogNormal(beta, xi) (Jacobian included)", {
  out <- lkj_prior_chain(n_iter = 20000, seed = 11)
  keep <- seq(2001, 20000)

  log_sigma_occ <- log(out$sigma[1, 1, keep]) # occupied cluster: MH step
  # Batch-means standard errors at these settings are ~0.05 (mean) and
  # ~0.03 (sd); the old kernel targets LogNormal(beta + xi^2, xi), whose
  # mean of log sigma sits ~1 above beta, far outside this bound.
  expect_lt(abs(mean(log_sigma_occ) - out$beta), 0.2)
  expect_lt(abs(sd(log_sigma_occ) - out$xi), 0.15)
})

test_that("R Metropolis step targets LKJ(eta) (P = 2: r ~ 2 Beta(eta, eta) - 1)", {
  eta <- 2
  out <- lkj_prior_chain(n_iter = 30000, eta = eta, seed = 12)
  r <- out$r[3001:30000, 1]

  expect_lt(abs(mean(r)), 0.05)
  expect_lt(abs(var(r) - 1 / (2 * eta + 1)), 0.03)
})

test_that("batch scale Metropolis step targets S_loc + InvGamma(rho, theta)", {
  out <- lkj_prior_chain(n_iter = 30000, seed = 13)
  S_excess <- out$S[1, 1, 3001:30000] - 1

  # rho = 3, theta = 1 (defaults): median of InvGamma(3, 1) is 1 / qgamma(.5, 3)
  expect_lt(abs(median(S_excess) - 1 / qgamma(0.5, shape = 3, rate = 1)), 0.05)
})

test_that("r_log_det holds log det(R) straight after initialisation", {
  out <- lkj_prior_chain(K = 3, n_iter = 5, eta = 3, seed = 14)
  expect_equal(out$r_log_det_init, out$r_log_det_true, tolerance = 1e-10)
})

test_that("an empty cluster's forced prior draw is (mu, Sigma) ~ joint prior", {
  out <- lkj_prior_chain(K = 2, n_iter = 20000, seed = 15)
  P <- 2
  n <- dim(out$mu)[3]

  # Cluster index 2 (1-based) is empty. Under the joint prior
  # mu | Sigma ~ N(mu_0, Sigma / kappa), so
  # kappa * (mu - mu_0)' Sigma^-1 (mu - mu_0) ~ chi^2_P exactly, and
  # successive sweeps are independent draws. Drawing mu from the
  # PREVIOUS Sigma and then redrawing Sigma independently breaks this.
  stat <- vapply(seq_len(n), function(i) {
    d <- out$mu[, 2, i] - out$mu_0
    Sigma <- out$cov[, (P + 1):(2 * P), i]
    out$kappa * drop(t(d) %*% solve(Sigma, d))
  }, numeric(1))

  expect_gt(ks.test(stat, "pchisq", df = P)$p.value, 1e-3)

  log_sigma_empty <- log(out$sigma[1, 2, ])
  expect_gt(ks.test(log_sigma_empty, "pnorm", mean = out$beta, sd = out$xi)$p.value, 1e-3)
})

test_that("empty clusters are never counted as Metropolis acceptances", {
  out <- lkj_prior_chain(K = 2, n_iter = 200, seed = 16)
  expect_equal(out$mu_count[2], 0)
  expect_equal(out$r_count[2], 0)
  expect_equal(out$sigma_count[2], 0)
  expect_gt(out$mu_count[1], 0)
})

# ---- driver-level checks ----------------------------------------------------

# Labels 1, 2 only and all items fixed, with K_max = 3: cluster 3 is empty
# for every sweep.
empty_cluster_fit <- function(type, ...) {
  d <- make_characterization_data(N = 45, K = 2)
  set.seed(77)
  batchSemiSupervisedMixtureModel(
    d$X, n_iter = 60, thin = 5, d$labels, rep(1L, d$N), d$batch_vec,
    type = type, K_max = 3, verbose = FALSE,
    control = batchmixControl(auto_tune = FALSE), ...
  )
}

test_that("reported acceptance rates are zero for an always-empty cluster", {
  fit_mvn <- empty_cluster_fit("MVN")
  expect_equal(fit_mvn$cov_acceptance_rate[3], 0)
  expect_equal(fit_mvn$mu_acceptance_rate[3], 0)

  fit_mvt <- empty_cluster_fit("MVT")
  expect_equal(fit_mvt$cov_acceptance_rate[3], 0)
  expect_equal(fit_mvt$mu_acceptance_rate[3], 0)
  expect_equal(fit_mvt$t_df_acceptance_rate[3], 0)

  fit_lkj <- empty_cluster_fit("MVN_LKJ")
  expect_equal(fit_lkj$r_acceptance_rate[3], 0)
  expect_equal(fit_lkj$sigma_acceptance_rate[3], 0)
  expect_equal(fit_lkj$mu_acceptance_rate[3], 0)
})

# R-level windows are "larger = wider"; C++ receives their reciprocals
# (larger = tighter) and adapts those. A window that is far too tight has
# acceptance ~ 1, so adaptation must LOOSEN it, i.e. reduce the C++ value.
tuned_fit <- function(type, control_args) {
  d <- make_characterization_data(N = 45, K = 2)
  set.seed(78)
  batchSemiSupervisedMixtureModel(
    d$X, n_iter = 300, thin = 10, d$labels, d$fixed_none, d$batch_vec,
    type = type, K_max = 2, verbose = FALSE,
    control = do.call(batchmixControl, c(control_args, list(auto_tune = TRUE, n_burn = 250)))
  )
}

test_that("auto-tuning loosens a far-too-tight S window and tightens a far-too-wide one", {
  tight <- tuned_fit("MVN", list(S_proposal_window = 1e-5)) # C++ value 1e5
  expect_lt(tight$final_S_proposal_window, 1e5)

  wide <- tuned_fit("MVN", list(S_proposal_window = 100)) # C++ value 0.01
  expect_gt(wide$final_S_proposal_window, 0.01)
})

test_that("auto-tuning loosens a far-too-tight covariance window (MVN, MVT)", {
  for (type in c("MVN", "MVT")) {
    fit <- tuned_fit(type, list(cov_proposal_window = 1e-4)) # Wishart df 1e4
    expect_lt(fit$final_cov_proposal_window, 1e4)
  }
})

test_that("auto-tuning loosens a far-too-tight t degrees-of-freedom window (MVT)", {
  fit <- tuned_fit("MVT", list(t_df_proposal_window = 1e-4))
  expect_lt(fit$final_t_df_proposal_window, 1e4)
})

test_that("the covariance window never adapts below the Wishart's minimum df (= P)", {
  fit <- tuned_fit("MVN", list(cov_proposal_window = 1 / 3)) # df 3, P = 2
  expect_gte(fit$final_cov_proposal_window, 2)
})
