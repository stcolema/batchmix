# Partial pooling (of the batch weights and of the batch-scale concentration)
# is the default. These tests pin the default resolution, that the old
# behaviour is still one argument away, the on-load notice, and that the
# prior predictive simulator follows the same default model.

fit_default <- function(B = 2, ...) {
  d <- make_characterization_data(N = 60, K = 2, B = B)
  set.seed(101)
  runBatchMix(
    d$X, n_iter = 200, thin = 5, batch_vec = d$batch_vec, type = "MVN",
    K_max = 2, initial_labels = d$labels, fixed = d$fixed_none, verbose = FALSE,
    control = batchmixControl(n_burn = 100), ...
  )
}

test_that("defaults are partial pooling of the weights and estimation of the scale concentration", {
  fit <- fit_default(B = 3)
  expect_identical(fit$batch_weight_prior, "partial_pooling")
  expect_true(fit$sample_s_scale)

  # rho is genuinely estimated: its trace moves, and stays above 2.
  rho <- as.vector(fit$rho)
  expect_gt(length(unique(rho)), 1)
  expect_true(all(rho > 2))

  # Batches have their own weights, which differ within a draw.
  w_batch <- fit$w_batch
  expect_equal(dim(w_batch)[1:2], c(3, 2))
  expect_true(any(abs(w_batch[1, 1, ] - w_batch[2, 1, ]) > 1e-8))
})

test_that("a single batch falls back to global weights", {
  d <- make_characterization_data(N = 40, K = 2, B = 1)
  set.seed(102)
  fit <- runBatchMix(
    d$X, n_iter = 40, thin = 5, batch_vec = rep(1L, d$N), type = "MVN",
    K_max = 2, initial_labels = d$labels, fixed = d$fixed_none, verbose = FALSE
  )
  expect_identical(fit$batch_weight_prior, "global")
})

test_that("the previous behaviour is recovered explicitly", {
  fit <- fit_default(B = 3, batch_weight_prior = "global", sample_s_scale = FALSE)
  expect_identical(fit$batch_weight_prior, "global")
  expect_false(fit$sample_s_scale)
  expect_equal(length(unique(as.vector(fit$rho))), 1)
  # Every batch shares one weight vector.
  expect_equal(fit$w_batch[1, , 1], fit$w_batch[2, , 1])
})

test_that("fitBatchMix() forwards sample_s_scale and batch_weight_prior", {
  d <- make_characterization_data(N = 40, K = 2, B = 2)
  set.seed(103)
  chains <- suppressWarnings(fitBatchMix(
    d$X, n_chains = 2, n_iter = 40, thin = 5, batch_vec = d$batch_vec, type = "MVN",
    K_max = 2, verbose = FALSE, sample_s_scale = FALSE, batch_weight_prior = "global"
  ))
  expect_false(chains[[1]]$sample_s_scale)
  expect_identical(chains[[1]]$batch_weight_prior, "global")
})

test_that("an invalid batch_weight_prior is rejected", {
  d <- make_characterization_data(N = 30, K = 2, B = 2)
  expect_error(
    runBatchMix(d$X, n_iter = 10, thin = 5, batch_vec = d$batch_vec, type = "MVN",
      K_max = 2, initial_labels = d$labels, fixed = d$fixed_none, verbose = FALSE,
      batch_weight_prior = "nonsense")
  )
})

test_that("the on-load message warns of the breaking changes and points to the vignettes", {
  on_attach <- getFromNamespace(".onAttach", "batchmix")
  text <- paste(
    testthat::capture_messages(on_attach("", "batchmix")),
    collapse = " "
  )
  for (pattern in c(
    "partial_pooling", "sample_s_scale", "global", "batch_corrected_data",
    "calcBICM", "salso", "minVI", "browseVignettes", "batchmix_workflow",
    "batch_weight_priors", "covariance_models", "suppressPackageStartupMessages"
  )) {
    expect_match(text, pattern, fixed = TRUE, info = pattern)
  }
})

test_that("every vignette named by the on-load message exists", {
  vignette_dir <- testthat::test_path("..", "..", "vignettes")
  skip_if_not(dir.exists(vignette_dir))
  for (v in c("batchmix_workflow", "batch_weight_priors", "covariance_models")) {
    expect_true(file.exists(file.path(vignette_dir, paste0(v, ".Rmd"))), info = v)
  }
})

# ---- prior predictive simulator follows the default model -------------------

test_that("simulatePriorPredictive() draws batch-specific weights and an estimated rho by default", {
  set.seed(7)
  N <- 60; P <- 2; K <- 3
  X <- matrix(rnorm(N * P), N, P)
  batch_vec <- sample(0:2, N, replace = TRUE)

  sims <- simulatePriorPredictive(X, batch_vec, K = K, type = "MVN", n_datasets = 400)

  w_batch <- lapply(sims, function(s) s$params$w_batch)
  expect_equal(dim(w_batch[[1]]), c(3, K))
  expect_equal(rowSums(w_batch[[1]]), rep(1, 3), tolerance = 1e-10)
  expect_gt(max(abs(w_batch[[1]][1, ] - w_batch[[1]][2, ])), 1e-6)
  # Symmetric prior over clusters: the average weight is 1 / K.
  expect_lt(abs(mean(vapply(w_batch, function(w) w[1, 1], numeric(1))) - 1 / K), 0.07)

  # rho - 2 ~ Gamma(2, 1) a priori: mean rho = 4; the batch-scale prior mean
  # theta / (rho - 1) = 0.5 is held fixed whatever rho is drawn.
  rho <- vapply(sims, function(s) s$params$rho, numeric(1))
  expect_true(all(rho > 2))
  expect_lt(abs(mean(rho) - 4), 0.4)
  S_excess <- unlist(lapply(sims, function(s) s$params$S - 1))
  expect_true(all(S_excess > 0))
  expect_lt(abs(mean(S_excess) - 0.5), 0.15)
})

test_that("simulatePriorPredictive() reproduces the fixed-hyperparameter model on request", {
  set.seed(8)
  N <- 40; P <- 2; K <- 2
  X <- matrix(rnorm(N * P), N, P)
  batch_vec <- sample(0:1, N, replace = TRUE)

  sims <- simulatePriorPredictive(
    X, batch_vec, K = K, type = "MVN", n_datasets = 5,
    batch_weight_prior = "global", sample_s_scale = FALSE, m_scale = 0.01
  )
  for (s in sims) {
    expect_equal(s$params$w_batch[1, ], s$params$w_batch[2, ])
    expect_equal(s$params$rho, 3)
  }
})

# ---- properties of the partial-pooling weight prior that the docs describe ----

test_that("the diffuse default pp_mu_prior_sd empties surplus components; a sharp one does not", {
  # Documented empirical behaviour (inst/experiments/overfitted_weight_priors.R):
  # not a theorem, so this compares two settings with a wide margin rather than
  # pinning a number.
  set.seed(11)
  N <- 240
  X <- rbind(cbind(rnorm(N / 2, 0), rnorm(N / 2, 0)), cbind(rnorm(N / 2, 7), rnorm(N / 2, 7)))
  batch_vec <- sample(1:2, N, replace = TRUE)
  mean_occupied <- function(sd, seed) {
    set.seed(seed)
    out <- suppressWarnings(runBatchMix(
      X, n_iter = 2000, thin = 20, batch_vec = batch_vec, type = "MVN", K_max = 6,
      verbose = FALSE, batch_weight_prior = "partial_pooling", pp_mu_prior_sd = sd
    ))
    lab <- out$samples[(nrow(out$samples) / 2 + 1):nrow(out$samples), , drop = FALSE]
    mean(apply(lab, 1, function(l) length(unique(l))))
  }
  diffuse <- mean(vapply(1:2, function(s) mean_occupied(10, s), numeric(1)))
  sharp <- mean(vapply(1:2, function(s) mean_occupied(1, s), numeric(1)))
  expect_lt(diffuse + 1, sharp)
})

test_that("the weight prior is exchangeable over clusters: no class is a priori more likely to dominate", {
  set.seed(12)
  N <- 40; K <- 6
  X <- matrix(rnorm(N * 2), N, 2)
  batch_vec <- sample(0:1, N, replace = TRUE)
  sims <- simulatePriorPredictive(X, batch_vec, K = K, type = "MVN", n_datasets = 1500)
  largest <- unlist(lapply(sims, function(s) max.col(s$params$w_batch, ties.method = "first")))
  freq <- tabulate(largest, K) / length(largest)
  # Each class should be the largest about 1 / K = 0.167 of the time; the
  # earlier reference-class parameterisation gave the last class ~0.01.
  expect_true(all(abs(freq - 1 / K) < 0.04))
})
