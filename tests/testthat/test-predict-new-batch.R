# Tests for predictNewBatch() (R/predictNewBatch.R): the no-data prior/GP-
# predictive case (i) for every batch_weight_prior, and the data-conditioned
# composition-sampling case (ii) (X_new supplied) for every model type.

# A small, well-separated 2-class, multi-batch dataset with a known batch
# shift/coordinate structure - deliberately NOT the shared random fixtures
# in helper-fixtures.R, since those aren't separated enough to check
# predictive accuracy against a known ground truth.
make_predict_fixture <- function(seed, B = 4, P = 2, n_per_batch = 30,
                                  mu_true = rbind(c(0, 0), c(6, 6)), sd = 0.8,
                                  batch_coordinates = NULL) {
  set.seed(seed)
  N <- B * n_per_batch
  X <- matrix(0, N, P)
  labels_true <- integer(N)
  batch_vec <- integer(N)
  row_i <- 1
  for (b in seq_len(B)) {
    for (k in sample(0:1, n_per_batch, replace = TRUE)) {
      X[row_i, ] <- rnorm(P, mu_true[k + 1, ], sd)
      labels_true[row_i] <- k
      batch_vec[row_i] <- b - 1
      row_i <- row_i + 1
    }
  }
  init_means <- cbind(colMeans(X[labels_true == 0, , drop = FALSE]), colMeans(X[labels_true == 1, , drop = FALSE]))
  list(
    X = X, labels_true = labels_true, batch_vec = batch_vec, fixed = rep(0L, N),
    mu_true = mu_true, init_means = init_means, P = P, B = B,
    batch_coordinates = batch_coordinates %||% seq_len(B)
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

make_new_batch <- function(seed, mu_true, shift_true, n_new = 40, sd = 0.8) {
  set.seed(seed)
  labels_new_true <- sample(0:1, n_new, replace = TRUE)
  P <- ncol(mu_true)
  X_new <- matrix(0, n_new, P)
  for (i in seq_len(n_new)) {
    X_new[i, ] <- rnorm(P, mu_true[labels_new_true[i] + 1, ] + shift_true, sd)
  }
  list(X_new = X_new, labels_true = labels_new_true)
}

test_that("predictNewBatch(): no-data case gives valid weight/shift/scale draws for every batch_weight_prior", {
  d <- make_predict_fixture(seed = 5001)

  fit_and_predict <- function(prior, fit_args = list(), predict_args = list()) {
    out <- do.call(batchSemiSupervisedMixtureModel, c(list(
      d$X, n_iter = 100, thin = 5, d$labels_true, d$fixed, d$batch_vec, "MVN",
      initial_class_means = d$init_means, batch_weight_prior = prior, verbose = FALSE
    ), fit_args))
    proc <- processMCMCChain(out, burn = 50)
    do.call(predictNewBatch, c(list(proc, d$X), predict_args))
  }

  pred_global <- fit_and_predict("global")
  pred_pp <- fit_and_predict("partial_pooling")
  pred_gp <- fit_and_predict("gp",
    fit_args = list(batch_coordinates = d$batch_coordinates),
    predict_args = list(new_batch_coordinate = max(d$batch_coordinates) + 1)
  )

  for (pred in list(pred_global, pred_pp, pred_gp)) {
    expect_true(all(abs(rowSums(pred$weight_draws) - 1) < 1e-8))
    expect_true(all(pred$weight_draws >= 0 & pred$weight_draws <= 1))
    expect_true(all(is.finite(pred$shift_draws)))
    expect_true(all(pred$scale_draws > 1)) # S_loc = 1 is the scale's lower bound
    expect_equal(dim(pred$mean_sum_est), c(d$P, 2))
  }
})

test_that("predictNewBatch(): 'gp' requires new_batch_coordinate, others warn if given one anyway", {
  d <- make_predict_fixture(seed = 5002)
  out <- batchSemiSupervisedMixtureModel(
    d$X, n_iter = 60, thin = 5, d$labels_true, d$fixed, d$batch_vec, "MVN",
    initial_class_means = d$init_means, batch_weight_prior = "gp",
    batch_coordinates = d$batch_coordinates, verbose = FALSE
  )
  proc <- processMCMCChain(out, burn = 30)
  expect_error(predictNewBatch(proc, d$X), "new_batch_coordinate is required")

  out_global <- batchSemiSupervisedMixtureModel(
    d$X, n_iter = 60, thin = 5, d$labels_true, d$fixed, d$batch_vec, "MVN",
    initial_class_means = d$init_means, verbose = FALSE
  )
  proc_global <- processMCMCChain(out_global, burn = 30)
  expect_warning(predictNewBatch(proc_global, d$X, new_batch_coordinate = 1), "ignored")
})

test_that("predictNewBatch(): GP prediction has higher uncertainty extrapolating than interpolating", {
  d <- make_predict_fixture(seed = 5003, B = 6, n_per_batch = 25)
  out <- batchSemiSupervisedMixtureModel(
    d$X, n_iter = 300, thin = 5, d$labels_true, d$fixed, d$batch_vec, "MVN",
    initial_class_means = d$init_means, batch_weight_prior = "gp",
    batch_coordinates = d$batch_coordinates, sample_gp_hyperparameters = TRUE, verbose = FALSE
  )
  proc <- processMCMCChain(out, burn = 150)

  pred_near <- predictNewBatch(proc, d$X, new_batch_coordinate = mean(d$batch_coordinates))
  pred_far <- predictNewBatch(proc, d$X, new_batch_coordinate = max(d$batch_coordinates) + 20)

  expect_gt(var(pred_far$weight_draws[, 1]), var(pred_near$weight_draws[, 1]))
})

test_that("predictNewBatch(): X_new requires batch_vec, rejects include_interaction fits", {
  d <- make_predict_fixture(seed = 5004)
  out <- batchSemiSupervisedMixtureModel(
    d$X, n_iter = 60, thin = 5, d$labels_true, d$fixed, d$batch_vec, "MVN",
    initial_class_means = d$init_means, verbose = FALSE
  )
  proc <- processMCMCChain(out, burn = 30)
  X_new <- d$X[1:5, , drop = FALSE]
  expect_error(predictNewBatch(proc, d$X, X_new = X_new), "batch_vec")

  out_interact <- batchSemiSupervisedMixtureModel(
    d$X, n_iter = 60, thin = 5, d$labels_true, d$fixed, d$batch_vec, "MVN",
    initial_class_means = d$init_means, include_interaction = TRUE, verbose = FALSE
  )
  proc_interact <- processMCMCChain(out_interact, burn = 30)
  expect_error(
    predictNewBatch(proc_interact, d$X, batch_vec = d$batch_vec, X_new = X_new),
    "include_interaction"
  )
})

test_that("predictNewBatch(): composition sampling recovers an injected batch shift and classifies correctly (MVN)", {
  skip_on_cran()
  d <- make_predict_fixture(seed = 5101, B = 3, n_per_batch = 50)
  out <- batchSemiSupervisedMixtureModel(
    d$X, n_iter = 1500, thin = 5, d$labels_true, d$fixed, d$batch_vec, "MVN",
    initial_class_means = d$init_means, verbose = FALSE
  )
  proc <- processMCMCChain(out, burn = 750)

  shift_true <- c(-2, -2)
  nb <- make_new_batch(5102, d$mu_true, shift_true)

  pred <- predictNewBatch(proc, d$X, batch_vec = d$batch_vec, X_new = nb$X_new,
                           n_draws = 15, n_pred_iter = 150, pred_thin = 2)

  expect_true(all(pred$shift_est < 0)) # correct sign, given a well-separated -2,-2 injected shift
  expect_gt(mean((pred$pred - 1) == nb$labels_true), 0.85)
  expect_true(all(abs(rowSums(pred$allocation_probability) - 1) < 1e-6))
  expect_equal(nrow(pred$allocation_probability), length(nb$labels_true))
})

test_that(".predictiveShiftScaleDraws(): recomputes theta from rho when sample_s_scale = TRUE, pinning E[S_b - S_loc]", {
  # For S_b - S_loc ~ InvGamma(shape = rho, scale = theta) (mean = theta /
  # (rho - 1) for rho > 1), theta = s_scale_prior_mean * (rho - 1) makes
  # that mean EXACTLY s_scale_prior_mean for any rho - the reparametrisation
  # sScaleConcentrationMetropolis() relies on to free only the concentration,
  # never the prior mean (see its documentation in src/mvnSampler.cpp). This
  # holds even when rho varies draw-to-draw, which a naive
  # `rep(processed_chain$theta, n_saved)` (the pre-partial-pooling behaviour,
  # still correct when sample_s_scale = FALSE) would not reproduce.
  P <- 2
  n_saved <- 4000
  rho_fixed <- 5 # constant across draws here - isolates the mean/variance
  # check from any confound with the Monte Carlo error of drawing rho itself
  s_scale_prior_mean <- 2.0
  processed_chain <- list(
    P = P,
    means = array(0, dim = c(P, 2, n_saved)),
    m_scale = 0.01, sample_m_scale = FALSE,
    rho = rep(rho_fixed, n_saved), theta = 999, # theta must be ignored below
    sample_s_scale = TRUE, s_scale_prior_mean = s_scale_prior_mean
  )
  X <- matrix(stats::rnorm(20 * P), ncol = P)

  set.seed(20240601)
  out <- batchmix:::.predictiveShiftScaleDraws(processed_chain, X)
  expect_equal(dim(out$scale), c(P, n_saved))

  excess <- out$scale - 1.0 # S_loc = 1.0
  # Monte Carlo SE of the mean, for InvGamma(rho, theta) with
  # theta = s_scale_prior_mean * (rho - 1): Var = theta^2 / ((rho-1)^2 (rho-2))
  # = s_scale_prior_mean^2 / (rho - 2).
  theoretical_var <- s_scale_prior_mean^2 / (rho_fixed - 2)
  mc_se <- sqrt(theoretical_var / (P * n_saved))
  expect_equal(mean(excess), s_scale_prior_mean, tolerance = 6 * mc_se / s_scale_prior_mean)

  # sample_s_scale = FALSE: theta_draws must fall back to the fixed,
  # originally-fitted scalar theta (rho's own value is then irrelevant to
  # theta, exactly the pre-partial-pooling behaviour).
  processed_chain_off <- processed_chain
  processed_chain_off$sample_s_scale <- FALSE
  processed_chain_off$theta <- 3.0
  set.seed(20240601)
  out_off <- batchmix:::.predictiveShiftScaleDraws(processed_chain_off, X)
  expect_equal(mean(out_off$scale - 1.0), 3.0 / (rho_fixed - 1), tolerance = 6 * sqrt((3.0^2 / ((rho_fixed - 1)^2 * (rho_fixed - 2))) / (P * n_saved)) / (3.0 / (rho_fixed - 1)))
})

test_that("predictNewBatch(): composition sampling conditions the new batch's scale on rho/theta when sample_s_scale = TRUE (MVN)", {
  skip_on_cran()
  d <- make_predict_fixture(seed = 5301, B = 3, n_per_batch = 50)
  out <- batchSemiSupervisedMixtureModel(
    d$X, n_iter = 1500, thin = 5, d$labels_true, d$fixed, d$batch_vec, "MVN",
    initial_class_means = d$init_means, sample_s_scale = TRUE, verbose = FALSE
  )
  expect_true(is.numeric(out$rho))
  expect_gt(length(out$rho), 1) # a per-iteration trace, not the scalar constructor argument
  expect_true(out$sample_s_scale)

  proc <- processMCMCChain(out, burn = 750)
  expect_true(isTRUE(proc$sample_s_scale))
  expect_true(is.numeric(proc$s_scale_prior_mean) && length(proc$s_scale_prior_mean) == 1)

  shift_true <- c(-2, -2)
  nb <- make_new_batch(5302, d$mu_true, shift_true)

  pred <- predictNewBatch(proc, d$X, batch_vec = d$batch_vec, X_new = nb$X_new,
                           n_draws = 15, n_pred_iter = 150, pred_thin = 2)

  expect_true(all(is.finite(pred$scale_new_draws)))
  expect_true(all(pred$scale_new_draws > 1)) # S_loc = 1 is the scale's lower bound
  expect_gt(mean((pred$pred - 1) == nb$labels_true), 0.8)
})

test_that("predictNewBatch(): composition sampling works for MVT/MVN_LKJ/MVN_MIXED", {
  skip_on_cran()

  for (type in c("MVT", "MVN_LKJ", "MVN_MIXED")) {
    d <- make_predict_fixture(seed = 5201, B = 3, n_per_batch = 40)
    shift_true <- c(-2, -2)
    nb <- make_new_batch(5202, d$mu_true, shift_true)

    fit_args <- list(
      X = d$X, n_iter = 1200, thin = 5, initial_labels = d$labels_true, fixed = d$fixed,
      batch_vec = d$batch_vec, type = type, initial_class_means = d$init_means, verbose = FALSE
    )
    pred_args <- list(
      X = d$X, batch_vec = d$batch_vec, X_new = nb$X_new,
      n_draws = 12, n_pred_iter = 120, pred_thin = 2
    )
    if (type == "MVN_MIXED") {
      column_type <- rep(0L, d$P)
      censor_code <- matrix(0L, nrow(d$X), d$P)
      fit_args$column_type <- column_type
      fit_args$censor_code <- censor_code
      pred_args$column_type <- column_type
      pred_args$censor_code <- censor_code
    }

    out <- do.call(batchSemiSupervisedMixtureModel, fit_args)
    proc <- processMCMCChain(out, burn = 600)
    pred_args$processed_chain <- proc
    pred <- do.call(predictNewBatch, pred_args)

    expect_gt(mean((pred$pred - 1) == nb$labels_true), 0.8, label = paste(type, "accuracy"))
    expect_true(all(abs(rowSums(pred$allocation_probability) - 1) < 1e-6), label = paste(type, "alloc sums to 1"))
  }
})
