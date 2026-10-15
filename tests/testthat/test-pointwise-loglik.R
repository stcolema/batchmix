# pointwiseLogLik() (see R/pointwiseLogLik.R). The sampler records the
# observed-data log-likelihood of every draw, which must equal the row sums of the
# pointwise values; that is the main oracle here, alongside independent
# recomputations from the cluster and batch parameters.

pw_data <- function(N = 60, P = 2) {
  set.seed(3)
  X <- rbind(
    matrix(rnorm(N / 2 * P, 0, 1), ncol = P),
    matrix(rnorm(N / 2 * P, 3, 1), ncol = P)
  )
  list(X = X, batch_vec = rep(1:3, length.out = nrow(X)))
}

pw_fit <- function(d, type = "MVN", ...) {
  runBatchMix(d$X, n_iter = 100, thin = 5, batch_vec = d$batch_vec, type = type, K_max = 3, verbose = FALSE, ...)
}

# The sampler records the sum over items, retained from the first draw after burn in
recorded <- function(fit, burn = 50) as.numeric(fit$observed_likelihood)[seq(floor(burn / fit$thin) + 1, nrow(fit$samples))]

test_that("row sums equal the observed likelihood the sampler records, for each model variant", {
  d <- pw_data()
  variants <- list(
    mvn_partial_pooling = list(type = "MVN"),
    mvn_global = list(type = "MVN", batch_weight_prior = "global"),
    mvn_gp = list(type = "MVN", batch_weight_prior = "gp", batch_coordinates = c(0, 1, 2)),
    mvn_interaction = list(type = "MVN", include_interaction = TRUE),
    mvt = list(type = "MVT"),
    mvn_lkj = list(type = "MVN_LKJ")
  )
  for (nm in names(variants)) {
    fit <- do.call(pw_fit, c(list(d), variants[[nm]]))
    ll <- pointwiseLogLik(fit, d$X, d$batch_vec, burn = 50)
    expect_equal(rowSums(ll), recorded(fit), tolerance = 1e-8, info = nm)
  }
})

test_that("semi-supervised items contribute the joint density at their label", {
  d <- pw_data()
  N <- nrow(d$X)
  labels <- c(rep(1, 10), rep(2, 10), rep(1, N - 20))
  fixed <- c(rep(1, 20), rep(0, N - 20))
  fit <- runBatchMix(d$X,
    n_iter = 100, thin = 5, batch_vec = d$batch_vec, type = "MVN", K_max = 3,
    initial_labels = labels, fixed = fixed, verbose = FALSE
  )
  expect_error(pointwiseLogLik(fit, d$X, d$batch_vec), "semi-supervised")
  ll <- pointwiseLogLik(fit, d$X, d$batch_vec, burn = 50, fixed = fixed)
  expect_equal(rowSums(ll), recorded(fit), tolerance = 1e-8)
  # a labelled item lies below the marginal over labels, which sums every cluster
  marginal <- pointwiseLogLik(fit, d$X, d$batch_vec, burn = 50, fixed = rep(0, N))
  expect_true(all(ll[, 1:20] <= marginal[, 1:20] + 1e-12))
  expect_equal(ll[, 21:N], marginal[, 21:N])
})

test_that("the value matches an independent computation from the cluster and batch parameters", {
  d <- pw_data()
  fit <- pw_fit(d)
  ll <- pointwiseLogLik(fit, d$X, d$batch_vec, burn = 50)
  P <- 2
  K <- 3
  b0 <- as.numeric(as.factor(d$batch_vec)) - 1
  reference <- function(r) {
    vapply(seq_len(nrow(d$X)), function(n) {
      b <- b0[n]
      terms <- vapply(0:(K - 1), function(k) {
        mu <- fit$means[, k + 1, r] + fit$batch_shift[, b + 1, r]
        Sk <- matrix(fit$covariance[, (k * P + 1):(k * P + P), r], P)
        # the batch scale inflates the diagonal of the cluster covariance only
        S <- Sk
        diag(S) <- diag(Sk) + (fit$batch_scale[, b + 1, r] - 1) * diag(Sk)
        dev <- d$X[n, ] - mu
        log(fit$w_batch[b + 1, k + 1, r]) -
          0.5 * (P * log(2 * pi) + log(det(S)) + drop(t(dev) %*% solve(S) %*% dev))
      }, numeric(1))
      max(terms) + log(sum(exp(terms - max(terms))))
    }, numeric(1))
  }
  first <- floor(50 / fit$thin) + 1
  for (i in c(1, 8, nrow(ll))) {
    expect_equal(ll[i, ], reference(first + i - 1), tolerance = 1e-10)
  }
})

test_that("a univariate t mixture matches the t density", {
  set.seed(8)
  x <- matrix(c(rnorm(30, 0, 1), rt(30, 4) + 4), ncol = 1)
  batch_vec <- rep(1:2, length.out = 60)
  fit <- runBatchMix(x, n_iter = 100, thin = 5, batch_vec = batch_vec, type = "MVT", K_max = 2, verbose = FALSE)
  ll <- pointwiseLogLik(fit, x, batch_vec, burn = 50)
  b0 <- batch_vec - 1
  r <- floor(50 / 5) + 3
  i <- r - floor(50 / 5)
  ref <- vapply(seq_along(x), function(n) {
    b <- b0[n]
    terms <- vapply(0:1, function(k) {
      cell <- k * 2 + b
      sigma <- sqrt(fit$cov_comb[1, cell + 1, r])
      mu <- fit$mean_sum[1, cell + 1, r]
      log(fit$w_batch[b + 1, k + 1, r]) +
        stats::dt((x[n] - mu) / sigma, df = fit$t_df[r, k + 1], log = TRUE) - log(sigma)
    }, numeric(1))
    max(terms) + log(sum(exp(terms - max(terms))))
  }, numeric(1))
  expect_equal(ll[i, ], ref, tolerance = 1e-10)
})

test_that("missing entries are integrated out", {
  d <- pw_data()
  fit <- pw_fit(d)
  X_na <- d$X
  X_na[1, 2] <- NA
  X_na[2, ] <- NA
  ll <- pointwiseLogLik(fit, d$X, d$batch_vec, burn = 50)
  ll_na <- pointwiseLogLik(fit, X_na, d$batch_vec, burn = 50)
  # an item with no observed coordinate carries no information
  expect_true(all(ll_na[, 2] == 0))
  # the others are unchanged
  expect_equal(ll_na[, -(1:2)], ll[, -(1:2)])
  # item 1 keeps coordinate 1: the marginal normal density on it
  P <- 2
  K <- 3
  B <- 3
  b <- 0
  r <- floor(50 / 5) + 1
  terms <- vapply(0:(K - 1), function(k) {
    cell <- k * B + b
    log(fit$w_batch[b + 1, k + 1, r]) +
      stats::dnorm(X_na[1, 1], fit$mean_sum[1, cell + 1, r], sqrt(fit$cov_comb[1, cell * P + 1, r]), log = TRUE)
  }, numeric(1))
  expect_equal(ll_na[1, 1], max(terms) + log(sum(exp(terms - max(terms)))), tolerance = 1e-10)
})

test_that("the result has the shape loo expects, and burn and processed chains are handled", {
  d <- pw_data()
  rownames(d$X) <- paste0("item", seq_len(nrow(d$X)))
  fits <- fitBatchMix(d$X,
    n_chains = 3, n_iter = 100, thin = 5, batch_vec = d$batch_vec, type = "MVN",
    K_max = 3, verbose = FALSE
  )
  ll <- pointwiseLogLik(fits, d$X, d$batch_vec)
  # 20 saved draws, the first 10 dropped (burn = n_iter / 2), three chains
  expect_equal(dim(ll), c(30, 60))
  expect_equal(colnames(ll), rownames(d$X))
  expect_equal(attr(ll, "chain_id"), rep(1:3, each = 10))
  expect_equal(nrow(pointwiseLogLik(fits, d$X, d$batch_vec, burn = 80)), 12)
  expect_error(pointwiseLogLik(fits, d$X, d$batch_vec, burn = 100), "no saved samples")
  # a single chain
  expect_equal(dim(pointwiseLogLik(fits[[1]], d$X, d$batch_vec)), c(10, 60))
  # processing a chain applies the burn in and changes nothing else
  processed <- processMCMCChains(fits, burn = 50)
  expect_equal(unclass(pointwiseLogLik(processed, d$X, d$batch_vec)), unclass(ll), ignore_attr = TRUE)
})

test_that("unsupported or inconsistent input is refused", {
  d <- pw_data()
  fit <- pw_fit(d)
  expect_error(pointwiseLogLik(fit, d$X[, 1, drop = FALSE], d$batch_vec), "must be a numeric")
  expect_error(pointwiseLogLik(fit, d$X, d$batch_vec[-1]), "one entry per row")
  expect_error(pointwiseLogLik(fit, d$X, rep(1:2, length.out = nrow(d$X))), "batches")
  expect_error(pointwiseLogLik(list(a = 1), d$X, d$batch_vec), "must be the output")

  N <- 200
  set.seed(1)
  Z <- matrix(rnorm(N * 3), N, 3)
  Z[, 3] <- (Z[, 3] > 0) * 1
  mixed <- runBatchMix(Z,
    n_iter = 50, thin = 5, batch_vec = rep(1:2, length.out = N), type = "MVN_MIXED", K_max = 2,
    column_type = c(0L, 0L, 1L), censor_code = matrix(0L, N, 3), verbose = FALSE
  )
  expect_error(pointwiseLogLik(mixed, Z, rep(1:2, length.out = N)), "MVN_MIXED")
})

test_that("the pointwise values feed loo and waic", {
  skip_if_not_installed("loo")
  d <- pw_data()
  fits <- fitBatchMix(d$X,
    n_chains = 2, n_iter = 200, thin = 5, batch_vec = d$batch_vec, type = "MVN",
    K_max = 3, verbose = FALSE
  )
  ll <- pointwiseLogLik(fits, d$X, d$batch_vec)
  r_eff <- suppressWarnings(loo::relative_eff(exp(ll), chain_id = attr(ll, "chain_id")))
  fit_loo <- suppressWarnings(loo::loo(ll, r_eff = r_eff))
  expect_s3_class(fit_loo, "loo")
  expect_equal(nrow(fit_loo$pointwise), 60)
  expect_true(is.finite(fit_loo$estimates["elpd_loo", "Estimate"]))
  # in-sample fit is above the leave-one-out estimate
  expect_gt(sum(log(colMeans(exp(ll)))), fit_loo$estimates["elpd_loo", "Estimate"])
  expect_s3_class(suppressWarnings(loo::waic(ll)), "waic")
})
