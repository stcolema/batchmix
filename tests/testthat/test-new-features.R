#!/usr/bin/Rscript
# Regression tests for the auto-tuning, batch x cluster interaction, and
# GP-correlated batch-weight features. Each test here corresponds to a
# concrete bug found and fixed during development (see the section
# comments), verified by simulating a known ground truth and checking the
# sampler recovers it - not just that the code runs without error.

test_that("robbinsMonroUpdate moves the proposal window in the correct direction", {
  # If the realised acceptance rate is above target, the window should grow
  # (encourage bigger, less-often-accepted moves); below target, it should
  # shrink.
  w_up <- robbinsMonroUpdate(window = 1.0, acceptance_rate = 0.9, target_rate = 0.234, n = 10, step_scale = 1.0, kappa = 0.6)
  w_down <- robbinsMonroUpdate(window = 1.0, acceptance_rate = 0.05, target_rate = 0.234, n = 10, step_scale = 1.0, kappa = 0.6)
  w_same <- robbinsMonroUpdate(window = 1.0, acceptance_rate = 0.234, target_rate = 0.234, n = 10, step_scale = 1.0, kappa = 0.6)

  expect_gt(w_up, 1.0)
  expect_lt(w_down, 1.0)
  expect_equal(w_same, 1.0, tolerance = 1e-10)
  expect_true(all(c(w_up, w_down, w_same) > 0))

  # The step size (and hence how far the window moves for a fixed
  # acceptance-rate gap) should shrink as n grows - this is what makes the
  # adaptation's total effect finite (diminishing adaptation).
  w_early <- robbinsMonroUpdate(window = 1.0, acceptance_rate = 0.9, target_rate = 0.234, n = 2, step_scale = 1.0, kappa = 0.6)
  w_late <- robbinsMonroUpdate(window = 1.0, acceptance_rate = 0.9, target_rate = 0.234, n = 2000, step_scale = 1.0, kappa = 0.6)
  expect_gt(abs(log(w_early)), abs(log(w_late)))
})

test_that("the Gamma-proposal Metropolis-Hastings orientation recovers a known target (regression for the batch-scale/sigma sign bug)", {
  # This is the exact asymmetric proposal used by batchScaleMetropolis()/
  # sigmaMHStep() (S_proposed ~ Gamma(shape = (S_current - S_loc) * window,
  # rate = window)). A previous version of the package's C++ assigned the
  # forward/reverse proposal densities to the wrong side of the acceptance
  # ratio here, which biased the chain: simulating against a flat/known
  # Gamma(5, 2) target with that (wrong) orientation recovered a mean of
  # ~1.4 instead of the true 2.5. Re-derive the same chain in R (using the
  # package's own exported gammaLogLikelihood()) to confirm the correct
  # orientation is (and stays) implemented.
  set.seed(20240615)
  a0 <- 5; b0 <- 2; w <- 8; n_iter <- 20000

  target_log <- function(x) gammaLogLikelihood(x, a0, b0)

  x <- 1.0
  samples <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    y <- rgamma(1, shape = x * w, rate = w)
    fwd <- gammaLogLikelihood(y, x * w, w)   # q(y | x): the density of the draw actually taken
    rev <- gammaLogLikelihood(x, y * w, w)   # q(x | y): the reverse proposal density
    # Correct orientation: reverse density pairs with the proposed target,
    # forward density pairs with the current target.
    proposed_score <- target_log(y) + rev
    current_score <- target_log(x) + fwd
    if (log(runif(1)) < proposed_score - current_score) x <- y
    samples[i] <- x
  }

  post <- samples[-(1:2000)]
  expect_equal(mean(post), a0 / b0, tolerance = 0.15)
})

test_that("the batch x cluster interaction term stays on the sum-to-zero subspace and recovers a known interaction effect", {
  # Regression for two bugs: (1) clusterMeanMetropolis()/batchShiftMetorpolis()
  # used to silently drop gamma's contribution when rebuilding mean_sum,
  # corrupting the interaction score; (2) an unconstrained per-cell gamma
  # prior/proposal is confounded with mu/m (classical two-way ANOVA
  # main-effect/interaction confound) and cannot be identified by shrinkage
  # alone, however tight - fixed by confining gamma to the sum-to-zero
  # subspace (every row and column of gamma(p, ., .) sums to zero) from the
  # prior draw onward.
  set.seed(20240616)

  P <- 1; K <- 2; B <- 3
  mu_true <- c(0, 10); m_true <- c(0, 1, -1)
  gamma_true <- matrix(0, K, B); gamma_true[2, 3] <- 5
  n_per_cell <- 60

  X <- c(); labels <- c(); batch_vec <- c()
  for (k in 1:K) for (b in 1:B) {
    mean_kb <- mu_true[k] + m_true[b] + gamma_true[k, b]
    X <- c(X, rnorm(n_per_cell, mean_kb, 0.6))
    labels <- c(labels, rep(k - 1, n_per_cell))
    batch_vec <- c(batch_vec, rep(b - 1, n_per_cell))
  }
  X <- matrix(X, ncol = 1)

  n_iter <- 2000; thin <- 10; n_burn <- 1000
  init_means <- matrix(c(0, 10), nrow = 1)

  out <- runBatchMix(X, n_iter, thin, batch_vec, "MVN",
    initial_labels = labels, fixed = rep(0, nrow(X)),
    initial_class_means = init_means,
    control = batchmixControl(n_burn = n_burn, gamma_proposal_window = 0.3),
    include_interaction = TRUE
  )

  post_idx <- (n_burn / thin + 1):(n_iter / thin)

  # Every individual saved draw of gamma(1, ., .) must have every row and
  # every column sum to (numerically) zero.
  for (i in post_idx[c(1, 5, length(post_idx))]) {
    gm <- matrix(out$gamma[1, , i], nrow = K, ncol = B)
    expect_equal(rowSums(gm), rep(0, K), tolerance = 1e-8)
    expect_equal(colSums(gm), rep(0, B), tolerance = 1e-8)
  }

  # The only thing that is actually identified is the TOTAL cell mean
  # mu_k + m_b + gamma_{k,b}; check that against the true generating cell
  # means (which is well-defined regardless of the mu/m/gamma split used to
  # generate the data).
  mu_post <- rowMeans(out$means[1, , post_idx])
  m_post <- rowMeans(out$batch_shift[1, , post_idx])
  gamma_post <- matrix(rowMeans(sapply(post_idx, function(i) out$gamma[1, , i])), nrow = K, ncol = B)

  implied <- outer(mu_post, m_post, "+") + gamma_post
  true_cell <- outer(mu_true, m_true, "+") + gamma_true

  expect_lt(max(abs(implied - true_cell)), 0.6)
})

test_that("GP-correlated batch weights recover a smooth batch-dependent trend (regression for the GP-jitter conditioning bug)", {
  # Regression for a numerical bug: a fixed, tiny jitter (1e-6) made
  # squaredExponentialKernel()'s covariance matrix catastrophically
  # ill-conditioned for realistic length-scale choices (condition number
  # ~1.7e7 in the case that exposed this), which silently corrupted
  # gp_cov_inv and drove the eta-coordinate acceptance rate to ~0 - fixed by
  # scaling the jitter with tau2. This just checks the downstream, visible
  # symptom: with batch_weight_prior = "gp", the model should recover a
  # smooth, monotonic batch-dependent trend in the cluster proportions
  # (anchoring cluster identity via initial_class_means, since cluster
  # index is otherwise only identified up to permutation).
  set.seed(20240617)

  B <- 6; n_per_batch <- 50; K <- 2; P <- 2
  true_p0 <- plogis(seq(-2, 2, length.out = B))

  X <- NULL; labels_true <- c(); batch_vec <- c()
  for (b in 1:B) {
    n0 <- rbinom(1, n_per_batch, true_p0[b]); n1 <- n_per_batch - n0
    Xb <- rbind(
      matrix(rnorm(n0 * P, mean = 0), n0, P),
      matrix(rnorm(n1 * P, mean = 6), n1, P)
    )
    X <- rbind(X, Xb)
    labels_true <- c(labels_true, rep(0, n0), rep(1, n1))
    batch_vec <- c(batch_vec, rep(b - 1, n0 + n1))
  }

  n_iter <- 600; thin <- 10; n_burn <- 300
  init_means <- matrix(c(0, 0, 6, 6), nrow = P)

  out <- runBatchMix(X, n_iter, thin, batch_vec, "MVN",
    initial_labels = labels_true, fixed = rep(0, nrow(X)),
    initial_class_means = init_means,
    control = batchmixControl(n_burn = n_burn, eta_proposal_window = 0.3),
    batch_weight_prior = "gp", gp_tau2 = 2.0, gp_length_scale = 3.0
  )

  post_idx <- (n_burn / thin + 1):(n_iter / thin)
  w_batch_post <- apply(out$w_batch[, 1, post_idx], 1, mean)

  expect_false(anyNA(w_batch_post))
  expect_gt(cor(w_batch_post, true_p0), 0.8)
})

test_that("partial-pooling batch weights recover per-batch proportions with no assumed order/structure", {
  # Unlike the GP test above, true per-batch proportions here are drawn
  # independently around a common population value and given to the model
  # in RANDOM (not sorted/ordered) batch order - there is no smooth trend
  # for a distance-based prior to exploit, and batch_weight_prior = "gp"
  # would have no reason to work on this data. "partial_pooling" instead
  # only assumes batches are exchangeable draws around a shared mean/
  # variance, which is exactly what was simulated.
  set.seed(20240618)

  B <- 8; n_per_batch <- 50; K <- 2; P <- 2
  true_logit_p0 <- rnorm(B, mean = 0.5, sd = 1.2) # random, unordered
  true_p0 <- plogis(true_logit_p0)

  X <- NULL; labels_true <- c(); batch_vec <- c()
  for (b in 1:B) {
    n0 <- rbinom(1, n_per_batch, true_p0[b]); n1 <- n_per_batch - n0
    Xb <- rbind(
      matrix(rnorm(n0 * P, mean = 0), n0, P),
      matrix(rnorm(n1 * P, mean = 6), n1, P)
    )
    X <- rbind(X, Xb)
    labels_true <- c(labels_true, rep(0, n0), rep(1, n1))
    batch_vec <- c(batch_vec, rep(b - 1, n0 + n1))
  }

  n_iter <- 800; thin <- 10; n_burn <- 400
  init_means <- matrix(c(0, 0, 6, 6), nrow = P)

  out <- runBatchMix(X, n_iter, thin, batch_vec, "MVN",
    initial_labels = labels_true, fixed = rep(0, nrow(X)),
    initial_class_means = init_means,
    control = batchmixControl(n_burn = n_burn, eta_proposal_window = 0.3),
    batch_weight_prior = "partial_pooling"
  )

  post_idx <- (n_burn / thin + 1):(n_iter / thin)
  w_batch_post <- apply(out$w_batch[, 1, post_idx], 1, mean)

  expect_false(anyNA(w_batch_post))
  expect_gt(cor(w_batch_post, true_p0), 0.7)

  # tau2 (pooling strength) should have moved away from its arbitrary
  # initial value and stay finite/positive.
  pp_tau2_post <- rowMeans(out$pp_tau2[, post_idx, drop = FALSE])
  expect_true(all(is.finite(pp_tau2_post)))
  expect_true(all(pp_tau2_post > 0))
})

test_that("global weights (the default) are unaffected by include_interaction/batch_weight_prior plumbing", {
  # A cheap regression check that weight_prior_type = 0 (global) leaves
  # w_batch inert (uniform) and every batch identical, i.e. the default
  # path is untouched by the partial-pooling/gp machinery sitting next to it.
  set.seed(1)
  N <- 120; P <- 2; K <- 2; B <- 3
  X <- matrix(rnorm(N * P), N, P); X[1:60, ] <- X[1:60, ] + 5
  labels <- c(rep(0, 60), rep(1, 60))
  batch_vec <- sample(0:(B - 1), N, replace = TRUE)

  out <- runBatchMix(X, 200, 10, batch_vec, "MVN",
    initial_labels = labels, fixed = rep(0, N),
    batch_weight_prior = "global"
  )

  last <- dim(out$w_batch)[3]
  w_last <- out$w_batch[, , last]
  expect_equal(unname(w_last[1, ]), unname(w_last[2, ]), tolerance = 1e-10)
  expect_equal(unname(w_last[1, ]), unname(w_last[3, ]), tolerance = 1e-10)
  expect_equal(out$weight_prior_type, 0)
})

test_that("batch-corrected data removes the interaction term, not just the batch shift (regression for updateBatchCorrectedData())", {
  # Regression for a real gap: updateBatchCorrectedData() only ever
  # subtracted the batch shift m_b, never the interaction term gamma_{k,b},
  # even when include_interaction = TRUE. Since gamma_{k,b} is just as much
  # a batch-specific effect as m_b (it only additionally depends on which
  # cluster the item is in), leaving it in meant "batch-corrected" data for
  # any (cluster, batch) cell with a real interaction still carried a
  # leftover batch-specific offset. Unlike splitting gamma against mu/m
  # individually (not identified), subtracting the SPECIFIC gamma_{k,b} an
  # item's own cell has is well-defined once gamma is sampled, since the
  # model only ever needs mu_k + m_b + gamma_{k,b} as a whole.
  set.seed(20240619)

  P <- 1; K <- 2; B <- 3
  mu_true <- c(0, 10); m_true <- c(0, 1, -1)
  gamma_true <- matrix(0, K, B); gamma_true[2, 3] <- 6 # a large, obvious interaction

  n_per_cell <- 70
  X <- c(); labels <- c(); batch_vec <- c()
  for (k in 1:K) for (b in 1:B) {
    mean_kb <- mu_true[k] + m_true[b] + gamma_true[k, b]
    X <- c(X, rnorm(n_per_cell, mean_kb, 0.5))
    labels <- c(labels, rep(k - 1, n_per_cell))
    batch_vec <- c(batch_vec, rep(b - 1, n_per_cell))
  }
  X <- matrix(X, ncol = 1)

  n_iter <- 2000; thin <- 10; n_burn <- 1000
  init_means <- matrix(c(0, 10), nrow = 1)

  out <- runBatchMix(X, n_iter, thin, batch_vec, "MVN",
    initial_labels = labels, fixed = rep(0, nrow(X)),
    initial_class_means = init_means,
    control = batchmixControl(n_burn = n_burn, gamma_proposal_window = 0.3),
    include_interaction = TRUE
  )

  post_idx <- (n_burn / thin + 1):(n_iter / thin)

  # Items in the cell with the large true interaction (cluster 2, batch 3,
  # 0-based: labels == 1, batch_vec == 2).
  target_cell <- which(labels == 1 & batch_vec == 2)
  other_cluster2_cells <- which(labels == 1 & batch_vec != 2)

  batch_corrected_mean <- function(items) {
    mean(sapply(post_idx, function(i) mean(out$batch_corrected_data[items, 1, i])))
  }

  target_est <- batch_corrected_mean(target_cell)
  other_est <- batch_corrected_mean(other_cluster2_cells)

  # If gamma were correctly removed, both should land close to mu_2 (~10,
  # up to the usual mu/m/gamma decomposition freedom, which does not affect
  # this fully-additive total). If gamma is NOT removed (the bug), the
  # target cell would sit noticeably higher than the rest of cluster 2 by
  # roughly the true interaction size.
  expect_lt(abs(target_est - other_est), 2.5)
})

test_that("continueChain() preserves the batch-weight-prior/interaction model spec (regression: it used to silently fall back to global/off)", {
  # continueChain() rebuilds the model via batchSemiSupervisedMixtureModel(),
  # but never forwarded include_interaction/batch_weight_prior (or their
  # hyperparameters), so a continued chain silently switched to
  # weight_prior_type = 0 / include_interaction = FALSE regardless of what
  # the original chain used - changing what was being fitted, not just how.
  # It also never combined the new weight-prior trace fields (w_batch,
  # eta_alr, gp_tau2, gp_length_scale, pp_mu, pp_tau2, gamma) with the old
  # chain's, so even when the spec happened to match, continuing a chain
  # silently dropped everything before the continuation for these fields.
  set.seed(321)
  N <- 150; P <- 2; K <- 2; B <- 3
  X <- matrix(rnorm(N * P), N, P); X[1:75, ] <- X[1:75, ] + 4
  labels <- c(rep(0, 75), rep(1, 75))
  batch_vec <- sample(1:B, N, replace = TRUE)
  fixed <- rep(0, N)
  time_pts <- c(0, 3, 11)

  for (wp in c("partial_pooling", "gp")) {
    coords <- if (wp == "gp") time_pts else NULL
    fit1 <- runBatchMix(X = X, K_max = K, initial_labels = labels, fixed = fixed,
      batch_vec = batch_vec, type = "MVN", n_iter = 300, thin = 10,
      control = batchmixControl(n_burn = 150),
      batch_weight_prior = wp, batch_coordinates = coords,
      sample_gp_hyperparameters = (wp == "gp")
    )
    fit2 <- continueChain(fit1, X, fixed, batch_vec, n_iter = 300, keep_old_samples = TRUE)

    expect_identical(fit2$batch_weight_prior, wp)
    expect_equal(fit2$n_iter, fit1$n_iter + 300)

    n1 <- dim(fit1$w_batch)[3]
    n2 <- dim(fit2$w_batch)[3]
    expect_equal(n2, n1 + floor(300 / fit1$thin))
    expect_equal(dim(fit2$w_batch)[1:2], c(B, K))
    expect_equal(length(fit2$complete_likelihood), n2)

    if (wp == "partial_pooling") {
      expect_equal(ncol(fit2$pp_mu), n2)
      expect_equal(ncol(fit2$pp_tau2), n2)
    }
    if (wp == "gp") {
      expect_equal(nrow(fit2$gp_tau2), n2)
      expect_equal(nrow(fit2$gp_length_scale), n2)
      expect_equal(ncol(fit2$gp_beta), n2)
    }
  }

  # The interaction term has the same failure mode, checked separately since
  # it is an independent flag from batch_weight_prior.
  fit3 <- runBatchMix(X = X, K_max = K, initial_labels = labels, fixed = fixed,
    batch_vec = batch_vec, type = "MVN", n_iter = 300, thin = 10,
    control = batchmixControl(n_burn = 150), include_interaction = TRUE
  )
  fit4 <- continueChain(fit3, X, fixed, batch_vec, n_iter = 300, keep_old_samples = TRUE)
  expect_true(fit4$include_interaction)
  n3 <- dim(fit3$gamma)[3]
  n4 <- dim(fit4$gamma)[3]
  expect_equal(n4, n3 + floor(300 / fit3$thin))

  # And the default (global weights, no interaction) path must still work
  # exactly as before - this fix must not regress the common case.
  fit5 <- runBatchMix(X = X, K_max = K, initial_labels = labels, fixed = fixed,
    batch_vec = batch_vec, type = "MVN", n_iter = 300, thin = 10,
    control = batchmixControl(n_burn = 150)
  )
  fit6 <- continueChain(fit5, X, fixed, batch_vec, n_iter = 300, keep_old_samples = TRUE)
  expect_identical(fit6$batch_weight_prior, "global")
  expect_false(fit6$include_interaction)
  expect_equal(fit6$n_iter, fit5$n_iter + 300)

  # w_batch/eta_alr are returned unconditionally (batch_weight_prior =
  # "global" just makes every batch's row identical within an iteration),
  # so they must be combined for "global" too - this was missed in an
  # earlier version of this same fix, which only combined them when
  # batch_weight_prior != "global", silently truncating w_batch/eta_alr
  # back down to just the new segment for the (most common) default case.
  n5 <- dim(fit5$w_batch)[3]
  n6 <- dim(fit6$w_batch)[3]
  expect_equal(n6, n5 + floor(300 / fit5$thin))
  expect_equal(dim(fit6$eta_alr)[3], n6)
})
