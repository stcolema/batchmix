# Non-stationary Gaussian-process batch-weight priors (gp_kernel = "rw1"/"rw2").

test_that("gpKernelMatrix: formulas, positive definiteness and non-stationarity", {
  x <- c(0, 0.2, 0.5, 1)
  K1 <- gpKernelMatrix(x, 1L, tau2 = 2, length_scale = 1, jitter = 0, level_var = 3, slope_var = 5)
  expect_equal(K1, 3 + 2 * outer(x, x, pmin))
  K2 <- gpKernelMatrix(x, 2L, tau2 = 2, length_scale = 1, jitter = 0, level_var = 3, slope_var = 5)
  m <- outer(x, x, pmin); M <- outer(x, x, pmax)
  expect_equal(K2, 3 + 5 * outer(x, x) + 2 * m^2 * (3 * M - m) / 6)
  # Variance grows with the coordinate (not stationary), unlike Matern-3/2.
  expect_true(all(diff(diag(K1)) > 0))
  expect_true(all(diff(diag(K2)) > 0))
  Km <- gpKernelMatrix(x, 0L, 2, 1, 0, 3, 5)
  expect_equal(diag(Km), rep(2, length(x)))
  for (K in list(K1, K2, Km)) {
    expect_true(min(eigen(K, symmetric = TRUE, only.values = TRUE)$values) > -1e-10)
  }
  expect_error(gpKernelMatrix(c(-1, 0), 1L, 1, 1, 0, 1, 1), "non-negative")
  expect_error(gpKernelMatrix(x, 3L, 1, 1, 0, 1, 1), "kernel type")
})

test_that("integrated Wiener prior has variance t^3/3 beyond its diffuse terms", {
  t <- c(0.3, 0.7, 1.0)
  K <- gpKernelMatrix(t, 2L, tau2 = 1, length_scale = 1, jitter = 0, level_var = 0, slope_var = 0)
  expect_equal(diag(K), t^3 / 3)
})

test_that("GP weight update samples each non-stationary prior (no data)", {
  x <- c(0, 0.2, 0.5, 1)
  for (kt in 1:2) {
    set.seed(100 + kt)
    out <- gpWeightChain(3L, 4L, 60000L, x, kt, tau2 = 1, length_scale = 1, level_sd = 1, eta_pw = 0.5)
    eta <- out$eta[, , 5001:60000]
    d <- eta[, 1, ] - eta[, 2, ] # level-free logit difference of two classes: covariance 2 Sigma
    target <- 2 * gpKernelMatrix(x, kt, 1, 1, 0, 1, 1)
    expect_lt(max(abs(cov(t(d)) - target)), 0.25)
    expect_lt(max(abs(rowMeans(d))), 0.1)
    rate <- out$eta_count[1] / 60000
    expect_gt(rate, 0.1)
    expect_lt(rate, 0.95)
  }
})

test_that("gp_kernel defaults to the non-stationary rw1 and stores the coordinate map", {
  set.seed(5)
  N <- 120; B <- 4
  batch_vec <- rep(1:B, each = N / B)
  labels <- sample(1:2, N, replace = TRUE)
  X <- matrix(rnorm(N * 2, mean = labels * 2.5), ncol = 2)
  coords <- c(2020, 2020.5, 2021, 2022)
  fit <- batchSemiSupervisedMixtureModel(
    X, n_iter = 60, thin = 1, labels, rep(0L, N), batch_vec, type = "MVN", K_max = 2,
    verbose = FALSE, batch_weight_prior = "gp", batch_coordinates = coords
  )
  expect_equal(fit$gp_kernel, "rw1")
  expect_equal(unique(as.numeric(fit$gp_tau2)), 20 / 3)
  expect_equal(fit$gp_coordinate_origin, 2020)
  expect_equal(fit$gp_coordinate_scale, 2)
  expect_equal(fit$batch_coordinates, coords)
  expect_true(all(is.finite(fit$w_batch)))

  fit_m <- batchSemiSupervisedMixtureModel(
    X, n_iter = 40, thin = 1, labels, rep(0L, N), batch_vec, type = "MVN", K_max = 2,
    verbose = FALSE, batch_weight_prior = "gp", batch_coordinates = coords, gp_kernel = "matern32"
  )
  expect_equal(fit_m$gp_kernel, "matern32")
  expect_equal(fit_m$gp_coordinate_scale, 1)
})

test_that("prediction widens with distance and rejects times before the origin", {
  set.seed(6)
  N <- 160; B <- 4
  batch_vec <- rep(1:B, each = N / B)
  labels <- sample(1:2, N, replace = TRUE)
  X <- matrix(rnorm(N * 2, mean = labels * 3), ncol = 2)
  coords <- c(0, 1, 2, 3)
  fit <- batchSemiSupervisedMixtureModel(
    X, n_iter = 200, thin = 1, labels, rep(0L, N), batch_vec, type = "MVN", K_max = 2,
    verbose = FALSE, batch_weight_prior = "gp", batch_coordinates = coords
  )
  pc <- processMCMCChain(fit, burn = 100)
  near <- batchmix:::.predictiveWeightDraws(pc, new_batch_coordinate = 3.2)$weight_draws
  far <- batchmix:::.predictiveWeightDraws(pc, new_batch_coordinate = 9)$weight_draws
  expect_gt(sd(far[, 1]), sd(near[, 1]))
  expect_error(batchmix:::.predictiveWeightDraws(pc, new_batch_coordinate = -1), "precedes")
})

test_that("predictNewBatch with data runs through the C++ composition sampler under rw2", {
  set.seed(8)
  N <- 160; B <- 4
  batch_vec <- rep(1:B, each = N / B)
  labels <- sample(1:2, N, replace = TRUE)
  X <- matrix(rnorm(N * 2, mean = labels * 3), ncol = 2)
  fit <- batchSemiSupervisedMixtureModel(
    X, n_iter = 120, thin = 1, labels, rep(0L, N), batch_vec, type = "MVN", K_max = 2,
    verbose = FALSE, batch_weight_prior = "gp", batch_coordinates = 0:3, gp_kernel = "rw2"
  )
  expect_equal(unique(as.numeric(fit$gp_tau2)), 500) # kernel-specific default start
  pc <- processMCMCChain(fit, burn = 60)
  lab_new <- sample(1:2, 40, replace = TRUE)
  X_new <- matrix(rnorm(80, mean = lab_new * 3), ncol = 2)
  pred <- predictNewBatch(pc, X = X, X_new = X_new, batch_vec = batch_vec,
    new_batch_coordinate = 4, n_draws = 5, n_pred_iter = 20)
  expect_true(all(is.finite(pred$weight_draws)))
})
