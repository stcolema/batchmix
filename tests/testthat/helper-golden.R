# Characterization ("golden master") test helper.
#
# These tests pin the CURRENT numeric output of batchSemiSupervisedMixtureModel()
# across all four model families (MVN, MVT, MVN_LKJ separation-strategy,
# MVN_MIXED) so that the planned collapse of the diamond-inherited
# mvnPredictive/mvtPredictive/mvnPredictiveSeparationStrategy/
# mvnPredictiveMixed/semisupervisedSampler classes can be checked for
# behavioural equivalence. Golden files live in tests/testthat/fixtures and
# are tracked in git deliberately (unlike testthat's own _snaps/, which is
# gitignored) - they are the empirical standard the refactor is checked
# against, not disposable snapshots.
#
# A mismatch here means the sampler's output for a fixed seed has changed.
# That is only expected immediately after a deliberate algorithmic change;
# when that happens the fixture must be regenerated deliberately (delete +
# re-run) and the numeric diff reviewed by hand before committing the new
# fixture.

expect_matches_golden <- function(value, name, dir = testthat::test_path("fixtures")) {
  path <- file.path(dir, paste0(name, ".rds"))

  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  newly_created <- !file.exists(path)
  if (newly_created) {
    saveRDS(value, path)
    testthat::skip(paste0(
      "Golden fixture '", name, "' did not exist; created at ", path,
      ". Re-run the tests to check against it, and commit the new file."
    ))
  }

  golden <- readRDS(path)
  testthat::expect_equal(value, golden, tolerance = 0)
}
