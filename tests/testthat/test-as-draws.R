# posterior::as_draws() methods (see R/asDraws.R)

draws_fit <- function(...) {
  set.seed(21)
  X <- rbind(matrix(rnorm(80, 0, 1), ncol = 2), matrix(rnorm(80, 3, 1), ncol = 2))
  batch_vec <- rep(1:3, length.out = nrow(X))
  list(
    X = X, batch_vec = batch_vec,
    fit = fitBatchMix(X,
      n_chains = 3, n_iter = 200, thin = 5, batch_vec = batch_vec,
      type = "MVN", K_max = 4, verbose = FALSE, ...
    )
  )
}

test_that("as_draws() gives a draws_array of the label-invariant quantities", {
  skip_if_not_installed("posterior")
  d <- draws_fit()
  dr <- posterior::as_draws(d$fit)
  expect_s3_class(dr, "draws_array")
  # 40 saved samples; the first 20 (burn = n_iter / 2 = 100, thin = 5) are dropped
  expect_equal(posterior::niterations(dr), 20)
  expect_equal(posterior::nchains(dr), 3)
  shifts <- as.vector(outer(1:2, 1:3, function(p, b) sprintf("batch_shift[%d,%d]", p, b)))
  scales <- as.vector(outer(1:2, 1:3, function(p, b) sprintf("batch_scale[%d,%d]", p, b)))
  expect_equal(
    posterior::variables(dr),
    c("log_lik_complete", "log_lik_observed", "bic", "occupied", "lambda_2", "rho", shifts, scales)
  )
  ch2 <- unclass(d$fit)[[2]]
  expect_equal(as.numeric(dr[, 2, "log_lik_complete"]), as.numeric(ch2$complete_likelihood)[21:40])
  expect_equal(as.numeric(dr[, 2, "bic"]), as.numeric(ch2$BIC)[21:40])
  expect_equal(as.numeric(dr[, 2, "rho"]), as.numeric(ch2$rho)[21:40])
  expect_equal(as.numeric(dr[, 2, "batch_shift[2,3]"]), ch2$batch_shift[2, 3, 21:40])
  expect_equal(as.numeric(dr[, 2, "batch_scale[1,2]"]), ch2$batch_scale[1, 2, 21:40])
  expect_equal(
    as.numeric(dr[, 2, "occupied"]),
    apply(ch2$samples[21:40, ], 1, function(z) length(unique(z)))
  )
})

test_that("the other formats and a single chain work, and burn is respected", {
  skip_if_not_installed("posterior")
  d <- draws_fit()
  expect_s3_class(posterior::as_draws_df(d$fit), "draws_df")
  expect_s3_class(posterior::as_draws_matrix(d$fit), "draws_matrix")
  expect_s3_class(posterior::as_draws_list(d$fit), "draws_list")
  expect_s3_class(posterior::as_draws_rvars(d$fit), "draws_rvars")
  expect_equal(nrow(posterior::as_draws_df(unclass(d$fit)[[1]])), 20)
  expect_equal(posterior::niterations(posterior::as_draws(d$fit, burn = 150)), 10)
  expect_error(posterior::as_draws(d$fit, burn = 200), "no saved samples")
})

test_that("processed and unprocessed chains give the same draws", {
  skip_if_not_installed("posterior")
  d <- draws_fit()
  raw <- posterior::as_draws(d$fit, burn = 100)
  pr <- posterior::as_draws(processMCMCChains(d$fit, burn = 100))
  expect_equal(unclass(raw), unclass(pr), ignore_attr = TRUE)
  # a processed chain with no burn in keeps every draw, as the default burn is not applied
  pr0 <- posterior::as_draws(processMCMCChains(d$fit, burn = 0))
  expect_equal(posterior::niterations(pr0), 40)
})

test_that("the posterior diagnostics agree with those of assessConvergence", {
  skip_if_not_installed("posterior")
  d <- draws_fit()
  s <- suppressWarnings(posterior::summarise_draws(
    posterior::as_draws(d$fit), "rhat", "ess_bulk", "ess_tail"
  ))
  a <- s[s$variable == "log_lik_complete", ]
  conv <- assessConvergence(d$fit, burn = 100)
  expect_equal(a$rhat, conv$rhat, tolerance = 1e-6)
  expect_equal(a$ess_bulk, conv$ess_bulk, tolerance = 1e-6)
  expect_equal(a$ess_tail, conv$ess_tail, tolerance = 1e-6)
})

test_that("sorted weights follow the batch weight prior, and allocations are added on request", {
  skip_if_not_installed("posterior")
  d <- draws_fit()
  dr <- posterior::as_draws(d$fit, weights = TRUE, allocations = TRUE, batch = FALSE)
  v <- posterior::variables(dr)
  expect_false(any(grepl("^batch_", v)))
  # partial pooling is the default: a set of weights for each of the 3 batches
  expect_true(all(as.vector(outer(1:3, 1:4, function(b, k) sprintf("sorted_weight[%d,%d]", b, k))) %in% v))
  expect_equal(sum(grepl("^allocation", v)), nrow(d$X))
  w <- posterior::subset_draws(dr, variable = sprintf("sorted_weight[2,%d]", 1:4))
  expect_true(all(abs(apply(w[, 1, ], 1, sum) - 1) < 1e-12))
  expect_true(all(apply(w[, 1, ], 1, function(r) !is.unsorted(rev(r)))))
  ch3 <- unclass(d$fit)[[3]]
  expect_equal(as.numeric(dr[, 3, "allocation[7]"]), ch3$samples[21:40, 7] + 1)
  # batch 2 of the third chain, not the pooled weights
  b2 <- t(apply(matrix(ch3$w_batch[2, , 21:40], ncol = 20), 2, function(r) sort(r / sum(r), decreasing = TRUE)))
  got <- matrix(as.numeric(dr[, 3, sprintf("sorted_weight[2,%d]", 1:4)]), nrow = 20)
  expect_equal(got, b2, ignore_attr = TRUE)

  g <- draws_fit(batch_weight_prior = "global")
  gv <- posterior::variables(posterior::as_draws(g$fit, weights = TRUE))
  expect_true(all(sprintf("sorted_weight[%d]", 1:4) %in% gv))
})

test_that("hyperparameters appear only when the sampler updated them", {
  skip_if_not_installed("posterior")
  d <- draws_fit(sample_s_scale = FALSE)
  expect_false("rho" %in% posterior::variables(posterior::as_draws(d$fit)))

  set.seed(4)
  X <- rbind(matrix(rnorm(60, 0, 1), ncol = 2), matrix(rnorm(60, 3, 1), ncol = 2))
  bv <- rep(1:4, length.out = nrow(X))
  gp_vars <- function(kernel) {
    fit <- fitBatchMix(X,
      n_chains = 2, n_iter = 100, thin = 5, batch_vec = bv, type = "MVN", K_max = 3,
      batch_weight_prior = "gp", batch_coordinates = 1:4, gp_kernel = kernel,
      sample_gp_hyperparameters = TRUE, verbose = FALSE
    )
    posterior::variables(posterior::as_draws(fit, batch = FALSE))
  }
  rw <- gp_vars("rw1")
  expect_true("gp_tau2" %in% rw)
  expect_false("gp_length_scale" %in% rw)
  expect_true(all(c("gp_tau2", "gp_length_scale") %in% gp_vars("matern32")))
})

test_that("chains with different settings are refused", {
  skip_if_not_installed("posterior")
  a <- unclass(draws_fit()$fit)[[1]]
  b <- unclass(draws_fit(sample_s_scale = FALSE)$fit)[[1]]
  expect_error(
    posterior::as_draws(structure(list(a, b), class = c("batchmix_fit_list", "list"))),
    "different settings"
  )
})
