# The batch-corrected data are the posterior mean of each item's latent
# batch-free signal: with x = mu_k + m_b + c + e, c ~ N(0, Sigma_k),
# e ~ N(0, diag((S_b - 1) * diag(Sigma_k))), Cov(x) = Sigma_kb (cov_comb) and
# E[c | x] = Sigma_k Sigma_kb^-1 (x - mu_k - m_b). See
# batchCorrectedPosteriorMean() in src/genericFunctions.cpp.

test_that("batchCorrectedPosteriorMean() is the exact posterior mean of the batch-free signal", {
  set.seed(1)
  P <- 2
  n <- 200000
  Sigma <- matrix(c(1, 0.6, 0.6, 1), P, P)
  S_b <- c(3, 2)
  mu <- c(1, -1)
  m_b <- c(0.5, 2)
  Sigma_kb <- Sigma + diag((S_b - 1) * diag(Sigma))

  c_true <- MASS::mvrnorm(n, rep(0, P), Sigma)
  e <- MASS::mvrnorm(n, rep(0, P), diag((S_b - 1) * diag(Sigma)))
  x <- sweep(c_true + e, 2, mu + m_b, "+")

  Y <- batchCorrectedPosteriorMean(
    X_t = t(x), labels = rep(0L, n), batch_vec = rep(0L, n), B = 1L,
    mu = matrix(mu, P, 1), mean_sum = matrix(mu + m_b, P, 1),
    cov = array(Sigma, c(P, P, 1)), cov_comb_inv = array(solve(Sigma_kb), c(P, P, 1))
  )
  c_hat <- sweep(Y, 2, mu, "-")

  # Defining properties of a conditional mean: the error is mean-zero and
  # orthogonal to the estimate.
  err <- c_true - c_hat
  expect_lt(max(abs(colMeans(err))), 0.02)
  expect_lt(max(abs(crossprod(err, c_hat) / n)), 0.02)

  # ...and it beats the residual/sqrt(S) rescaling on mean squared error,
  # which is neither unbiased nor orthogonal for a non-diagonal Sigma.
  c_rescaled <- sweep(sweep(x, 2, mu + m_b, "-"), 2, sqrt(S_b), "/")
  expect_lt(mean((c_true - c_hat)^2), mean((c_true - c_rescaled)^2))
  expect_gt(max(abs(crossprod(c_true - c_rescaled, c_rescaled) / n)), 0.05)
})

test_that("with no batch inflation (S = 1) the corrected data are the data less the shift", {
  P <- 3
  Sigma <- crossprod(matrix(c(1, 0.2, 0, 0.3, 1, 0.1, 0, 0.4, 1), P))
  x <- matrix(rnorm(5 * P), 5, P)
  mu <- c(0, 1, 2)
  m_b <- c(0.3, -0.2, 0.1)
  Y <- batchCorrectedPosteriorMean(
    t(x), rep(0L, 5), rep(0L, 5), 1L, matrix(mu, P, 1), matrix(mu + m_b, P, 1),
    array(Sigma, c(P, P, 1)), array(solve(Sigma), c(P, P, 1))
  )
  expect_equal(Y, sweep(x, 2, m_b, "-"), tolerance = 1e-10)
})

test_that("batch_corrected_data in every sampler's output matches the formula from the saved draws", {
  d <- make_characterization_data(N = 45, K = 2)

  for (type in c("MVN", "MVT", "MVN_LKJ")) {
    set.seed(31)
    fit <- batchSemiSupervisedMixtureModel(
      d$X, n_iter = 20, thin = 5, d$labels, d$fixed_none, d$batch_vec,
      type = type, K_max = 2, verbose = FALSE,
      control = batchmixControl(auto_tune = FALSE)
    )

    P <- d$P
    B <- d$B
    for (i in seq_len(dim(fit$batch_corrected_data)[3])) {
      expected <- t(vapply(seq_len(d$N), function(n) {
        k <- fit$samples[i, n] # 0-indexed cluster label at save time
        kb <- k * B + (d$batch_vec[n] - 1)
        cov_k <- fit$covariance[, (k * P + 1):((k + 1) * P), i]
        cov_kb <- fit$cov_comb[, (kb * P + 1):((kb + 1) * P), i]
        fit$means[, k + 1, i] +
          drop(cov_k %*% solve(cov_kb, fit$latent_data[n, , i] - fit$mean_sum[, kb + 1, i]))
      }, numeric(P)))
      expect_equal(fit$batch_corrected_data[, , i], expected, tolerance = 1e-8, ignore_attr = TRUE)
    }
  }
})

test_that("interaction term is removed from the corrected data along with the batch shift", {
  d <- make_characterization_data(N = 45, K = 2)
  set.seed(32)
  fit <- batchSemiSupervisedMixtureModel(
    d$X, n_iter = 20, thin = 5, d$labels, d$fixed_none, d$batch_vec,
    type = "MVN", K_max = 2, verbose = FALSE, include_interaction = TRUE,
    control = batchmixControl(auto_tune = FALSE)
  )

  P <- d$P
  B <- d$B
  i <- dim(fit$batch_corrected_data)[3]
  n <- 1
  k <- fit$samples[i, n]
  kb <- k * B + (d$batch_vec[n] - 1)
  cov_k <- fit$covariance[, (k * P + 1):((k + 1) * P), i]
  cov_kb <- fit$cov_comb[, (kb * P + 1):((kb + 1) * P), i]

  # mean_sum carries mu_k + m_b + gamma_{k,b}; the corrected value must
  # subtract all of it, not just mu_k + m_b.
  # gamma is saved flattened from a (P, K, B) cube, so column k + K * b.
  gamma_kb <- fit$gamma[, k + 2 * (d$batch_vec[n] - 1) + 1, i]
  expect_true(any(abs(gamma_kb) > 0))
  expect_equal(
    fit$mean_sum[, kb + 1, i],
    fit$means[, k + 1, i] + fit$batch_shift[, d$batch_vec[n], i] + gamma_kb,
    tolerance = 1e-8, ignore_attr = TRUE
  )
  expect_equal(
    drop(fit$batch_corrected_data[n, , i]),
    fit$means[, k + 1, i] +
      drop(cov_k %*% solve(cov_kb, fit$latent_data[n, , i] - fit$mean_sum[, kb + 1, i])),
    tolerance = 1e-8, ignore_attr = TRUE
  )
})
