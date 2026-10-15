# The namespace lists its exports explicitly; internal engines and helpers stay
# internal (they remain reachable to the other tests through the namespace).

test_that("the documented interface is exported", {
  exports <- getNamespaceExports("batchmix")
  expect_true(all(c(
    "fitBatchMix", "runBatchMix", "batchSemiSupervisedMixtureModel", "batchmixControl",
    "predictNewBatch", "gpKernelMatrix", "simulatePriorPredictive", "assessConvergence"
  ) %in% exports))
})

test_that("internal engines and helpers are not exported", {
  exports <- getNamespaceExports("batchmix")
  expect_false(any(c(
    "sampleSemisupervisedMVN", "predictNewBatchMVN", "gammaLogLikelihood",
    "wishartLogLikelihood", "robbinsMonroUpdate", "minVIGreedyRefine",
    "batchCorrectedPosteriorMean", "VI.lb"
  ) %in% exports))
})
