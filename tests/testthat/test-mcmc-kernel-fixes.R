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
    r_pw = 1.0, sigma_pw = 5, mu_pw = 2.5, m_pw = 0.3, S_pw = 5
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

# ---- efficiency changes: each must leave the target distribution unchanged ----

test_that("covariance-shaped mean random walk targets the prior (occupied cluster)", {
  out <- lkj_prior_chain(K = 2, n_iter = 30000, seed = 21)
  P <- 2
  keep <- 3001:30000
  stat <- vapply(keep, function(i) {
    d <- out$mu[, 1, i] - out$mu_0
    out$kappa * drop(t(d) %*% solve(out$cov[, 1:P, i], d))
  }, numeric(1))
  # chi^2_2: mean 2, median 2 log 2. Autocorrelated draws: loose bounds that a
  # mis-specified proposal correction (shifting the law) would break.
  expect_lt(abs(mean(stat) - P), 0.3)
  expect_lt(abs(median(stat) - qchisq(0.5, P)), 0.25)
})

test_that("interactionLogLikDelta() equals the full before/after log-likelihood difference", {
  set.seed(3)
  P <- 3; N <- 40; C <- 4
  X_t <- matrix(rnorm(P * N), P, N)
  cell <- sample(0:(C - 1), N, replace = TRUE)
  mean_sum <- matrix(rnorm(P * C), P, C)
  Lam <- array(0, c(P, P, C))
  for (c in seq_len(C)) Lam[, , c] <- solve(crossprod(matrix(rnorm(P * P), P)) + diag(P))
  p <- 1L
  delta <- rnorm(C, sd = 0.4)
  nu <- runif(N, 3, 20)

  loglik <- function(ms, t_df = NULL) {
    sum(vapply(seq_len(N), function(n) {
      c <- cell[n] + 1
      r <- X_t[, n] - ms[, c]
      u <- drop(t(r) %*% Lam[, , c] %*% r)
      if (is.null(t_df)) -0.5 * u else -0.5 * (t_df[n] + P) * log1p(u / t_df[n])
    }, numeric(1)))
  }
  shifted <- mean_sum
  shifted[p + 1, ] <- shifted[p + 1, ] + delta

  expect_equal(
    interactionLogLikDelta(X_t, cell, mean_sum, Lam, p, delta, numeric(0)),
    loglik(shifted) - loglik(mean_sum), tolerance = 1e-10
  )
  expect_equal(
    interactionLogLikDelta(X_t, cell, mean_sum, Lam, p, delta, nu),
    loglik(shifted, nu) - loglik(mean_sum, nu), tolerance = 1e-10
  )
})

test_that("partial-pooling update targets each batch's exact conditional (logit difference, K = 2)", {
  # Batch 0: 30 items, 20 in cluster 0; batch 1: 6 items, 1 in cluster 0;
  # batch 2: empty (prior only). For two classes the weights depend on the
  # logits only through d = eta_1 - eta_2, whose exact conditional given the
  # fixed hyperparameters is
  #   p(d) ~ exp(c_1 d - N log(1 + e^d)) N(d; mu_1 - mu_2, tau2_1 + tau2_2).
  labels <- c(rep(0L, 20), rep(1L, 10), 0L, rep(1L, 5))
  batch_vec <- c(rep(0L, 30), rep(1L, 6))
  mu <- c(0.4, -0.4); tau2 <- c(1.0, 0.8)

  set.seed(41)
  out <- partialPoolingWeightChain(labels, batch_vec, K = 2L, B = 3L, n_iter = 40000L,
    mu = mu, tau2 = tau2, eta_pw = 1.0)
  expect_equal(out$eta_moves_per_sweep, 3)

  d_chain <- out$eta[, 1, 5001:40000] - out$eta[, 2, 5001:40000] # B x draws
  quantiles_of_target <- function(c_b, n_b, probs) {
    grid <- seq(-12, 12, length.out = 40001)
    dens <- exp(c_b * grid - n_b * log1p(exp(grid)) +
      dnorm(grid, mu[1] - mu[2], sqrt(sum(tau2)), log = TRUE))
    cdf <- cumsum(dens) / sum(dens)
    vapply(probs, function(q) grid[which(cdf >= q)[1]], numeric(1))
  }
  probs <- c(0.1, 0.5, 0.9)
  for (spec in list(list(b = 1, c = 20, n = 30), list(b = 2, c = 1, n = 6), list(b = 3, c = 0, n = 0))) {
    expect_lt(
      max(abs(unname(quantile(d_chain[spec$b, ], probs)) - quantiles_of_target(spec$c, spec$n, probs))),
      0.15
    )
  }

  # Acceptance is a per-move rate in (0, 1), not a per-sweep count up to B.
  rate <- out$eta_count[1] / (40000 * out$eta_moves_per_sweep)
  expect_gt(rate, 0.05)
  expect_lt(rate, 0.95)
})

test_that("an empty batch's logits follow their prior exactly, for K = 4 classes", {
  # With no data the conditional is the prior: eta_bk ~ N(mu_k, tau2_k)
  # independently (the common shift is redrawn from its own exact Gibbs
  # conditional, so the shift-invariant logit differences are what is
  # identified).
  labels <- rep(0:3, each = 5)
  batch_vec <- rep(0L, 20) # batches 1 and 2 are empty
  mu <- c(1.0, 0.0, -0.5, -0.5); tau2 <- c(0.5, 1.0, 1.5, 2.0)

  set.seed(42)
  out <- partialPoolingWeightChain(labels, batch_vec, K = 4L, B = 3L, n_iter = 30000L,
    mu = mu, tau2 = tau2, eta_pw = 1.0)
  eta <- out$eta[3, , 3001:30000] # K x draws, an empty batch
  diffs <- eta[1, ] - eta[2, ]
  expect_lt(abs(mean(diffs) - (mu[1] - mu[2])), 0.06)
  expect_lt(abs(var(diffs) - (tau2[1] + tau2[2])), 0.15)
  diffs34 <- eta[3, ] - eta[4, ]
  expect_lt(abs(mean(diffs34)), 0.08)
  expect_lt(abs(var(diffs34) - (tau2[3] + tau2[4])), 0.25)
})

test_that("sample_s_scale = TRUE rejects rho <= 2 instead of silently freezing rho", {
  d <- make_characterization_data(N = 30, K = 2)
  expect_error(
    batchSemiSupervisedMixtureModel(
      d$X, n_iter = 10, thin = 5, d$labels, d$fixed_none, d$batch_vec,
      type = "MVN", K_max = 2, verbose = FALSE, rho = 1.5, sample_s_scale = TRUE
    ),
    "rho > 2"
  )
})
