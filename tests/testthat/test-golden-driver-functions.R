# Golden-master tests at the R-facing API (batchSemiSupervisedMixtureModel()),
# covering every model family that currently sits behind the diamond-
# inherited class hierarchy: MVN and MVT (mvnSampler/mvtSampler +
# mvnPredictive/mvtPredictive), MVN_LKJ (mvnSamplerSeparationStrategy +
# mvnPredictiveSeparationStrategy), and MVN_MIXED (mvnSamplerMixed +
# mvnPredictiveMixed). Each scenario with fixed_semi exercises both branches
# of updateAllocation() (fixed and unfixed items) in one run; the MVN
# scenario is additionally run fully unsupervised (fixed_none) to pin the
# fixed = all-zero path end to end.

test_that("MVN: unsupervised (fixed = all zero)", {
  d <- make_characterization_data()

  set.seed(1001)
  out <- batchSemiSupervisedMixtureModel(
    d$X, n_iter = 40, thin = 5, d$labels, d$fixed_none, d$batch_vec, type = "MVN",
    K_max = d$K, verbose = FALSE
  )

  expect_identical(out$Semisupervised, FALSE)
  expect_matches_golden(out, "mvn-unsupervised")
})

test_that("MVN: semi-supervised", {
  d <- make_characterization_data()

  set.seed(1002)
  out <- batchSemiSupervisedMixtureModel(
    d$X, n_iter = 40, thin = 5, d$labels, d$fixed_semi, d$batch_vec, type = "MVN",
    K_max = d$K, verbose = FALSE
  )

  expect_identical(out$Semisupervised, TRUE)
  expect_matches_golden(out, "mvn-semisupervised")
})

test_that("MVT: semi-supervised", {
  d <- make_characterization_data()

  set.seed(2001)
  out <- batchSemiSupervisedMixtureModel(
    d$X, n_iter = 40, thin = 5, d$labels, d$fixed_semi, d$batch_vec, type = "MVT",
    K_max = d$K, verbose = FALSE
  )

  expect_matches_golden(out, "mvt-semisupervised")
})

test_that("MVN_LKJ (separation-strategy): semi-supervised", {
  d <- make_characterization_data()

  set.seed(3001)
  out <- batchSemiSupervisedMixtureModel(
    d$X, n_iter = 40, thin = 5, d$labels, d$fixed_semi, d$batch_vec, type = "MVN_LKJ",
    K_max = d$K, verbose = FALSE
  )

  expect_matches_golden(out, "mvn-lkj-semisupervised")
})

test_that("MVN_MIXED: semi-supervised", {
  d <- make_mixed_characterization_data()

  set.seed(4001)
  out <- batchSemiSupervisedMixtureModel(
    d$X, n_iter = 40, thin = 5, d$labels, d$fixed_semi, d$batch_vec, type = "MVN_MIXED",
    K_max = d$K, verbose = FALSE,
    column_type = d$column_type, censor_code = d$censor_code
  )

  expect_matches_golden(out, "mvn-mixed-semisupervised")
})
