# Parallel chains: opt-in, reproducible streams, capped under R CMD check.

pc_data <- function(N = 40) {
  set.seed(5)
  list(
    X = rbind(matrix(rnorm(N, 0, 1), ncol = 2), matrix(rnorm(N, 3, 1), ncol = 2)),
    batch_vec = rep(1:2, length.out = N)
  )
}

pc_fit <- function(d, n_chains = 2, ...) {
  fitBatchMix(d$X, n_chains = n_chains, n_iter = 40, thin = 5, batch_vec = d$batch_vec, type = "MVN",
              K_max = 3, verbose = FALSE, ...)
}

with_env <- function(name, value, code) {
  old <- Sys.getenv(name, unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv(name) else do.call(Sys.setenv, stats::setNames(list(old), name)))
  do.call(Sys.setenv, stats::setNames(list(value), name))
  force(code)
}

test_that("the number of cores defaults to one and is capped by chains, cores and the check limit", {
  old <- options(batchmix.cores = NULL)
  on.exit(options(old))
  with_env("_R_CHECK_LIMIT_CORES_", "", {
    expect_equal(batchmix:::.batchmixResolveCores(NULL, 4), 1L)
    expect_equal(batchmix:::.batchmixResolveCores(1, 4), 1L)
    # never more cores than chains
    expect_equal(batchmix:::.batchmixResolveCores(64, 1), 1L)
    expect_lte(batchmix:::.batchmixResolveCores(64, 8), 8L)
    options(batchmix.cores = 2)
    expect_lte(batchmix:::.batchmixResolveCores(NULL, 4), 2L)
  })
  # as data.table does under R CMD check: at most two cores
  with_env("_R_CHECK_LIMIT_CORES_", "TRUE", {
    expect_lte(batchmix:::.batchmixResolveCores(64, 8), 2L)
  })
  for (bad in list(0, -1, NA, "a", c(1, 2))) {
    expect_error(batchmix:::.batchmixResolveCores(bad, 4), "n_cores")
  }
})

test_that("chain streams are reproducible from set.seed and leave the user's generator alone", {
  set.seed(10)
  a <- batchmix:::.batchmixStreamSeeds(3)
  set.seed(10)
  b <- batchmix:::.batchmixStreamSeeds(3)
  expect_identical(a, b)
  expect_length(unique(lapply(a, paste, collapse = " ")), 3)
  expect_true(all(vapply(a, function(s) s[1] %% 100 == 7, logical(1)))) # L'Ecuyer-CMRG kind code

  # The generator kind is unchanged, and its stream continues as if one draw had been taken
  set.seed(11)
  kind <- RNGkind()
  batchmix:::.batchmixStreamSeeds(3)
  expect_identical(RNGkind(), kind)
  after <- runif(1)
  set.seed(11)
  sample.int(.Machine$integer.max, 1L)
  expect_equal(after, runif(1))
})

test_that("serial runs (the default) are unchanged by the cores machinery", {
  d <- pc_data()
  set.seed(12)
  a <- pc_fit(d)
  set.seed(12)
  b <- lapply(1:2, function(i) {
    runBatchMix(d$X, 40, 5, d$batch_vec, "MVN", K_max = 3, verbose = FALSE)
  })
  for (i in 1:2) {
    expect_identical(a[[i]]$samples, b[[i]]$samples)
    expect_identical(a[[i]]$weights, b[[i]]$weights)
    expect_equal(a[[i]]$Chain, i)
  }
})

test_that("parallel chains are reproducible, differ from each other, and do not depend on the cores or the worker type", {
  skip_on_cran()
  d <- pc_data()
  run <- function(type, n_cores = 2) {
    old <- options(batchmix.parallel_type = type)
    on.exit(options(old))
    set.seed(13)
    pc_fit(d, n_chains = 3, n_cores = n_cores)
  }
  if (.Platform$OS.type != "windows") {
    f1 <- run("FORK")
    f2 <- run("FORK")
    f3 <- run("FORK", n_cores = 3)
    expect_s3_class(f1, "batchmix_fit_list")
    expect_length(f1, 3)
    expect_equal(vapply(f1, function(ch) ch$Chain, numeric(1)), 1:3)
    for (i in 1:3) {
      expect_identical(f1[[i]]$samples, f2[[i]]$samples)
      # the number of cores does not change the result
      expect_identical(f1[[i]]$samples, f3[[i]]$samples)
      expect_identical(f1[[i]]$weights, f3[[i]]$weights)
    }
    expect_false(identical(f1[[1]]$samples, f1[[2]]$samples))
    # the convergence report is still attached
    expect_false(is.null(attr(f1, "convergence")))
  }

  # Sockets need the installed package in the workers
  skip_if(nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_", "")) && .Platform$OS.type == "windows")
  sock <- run("PSOCK")
  expect_length(sock, 3)
  if (.Platform$OS.type != "windows") {
    for (i in 1:3) expect_identical(sock[[i]]$samples, f1[[i]]$samples)
  }
})

test_that("a failing chain is reported, and n_chains is checked", {
  skip_on_cran()
  skip_on_os("windows")
  d <- pc_data()
  expect_error(
    fitBatchMix(d$X, n_chains = 2, n_iter = 40, thin = 5, batch_vec = d$batch_vec[-1], type = "MVN",
                K_max = 3, n_cores = 2, verbose = FALSE),
    "failed"
  )
  for (bad in list(0, -1, NA, c(2, 3))) {
    expect_error(pc_fit(d, n_chains = bad), "n_chains")
  }
})

test_that("the cores are announced when verbose, and fitBatchMix passes n_cores on", {
  skip_on_cran()
  skip_on_os("windows")
  d <- pc_data()
  expect_message(
    fitBatchMix(d$X, n_chains = 2, n_iter = 40, thin = 5, batch_vec = d$batch_vec, type = "MVN",
                K_max = 3, n_cores = 2, verbose = TRUE),
    "on 2 cores"
  )
})

test_that("continueChains can continue the chains in parallel, reproducibly", {
  skip_on_cran()
  skip_on_os("windows")
  d <- pc_data()
  set.seed(14)
  fits <- pc_fit(d, n_chains = 2)
  cont <- function(n_cores) {
    set.seed(15)
    continueChains(fits, d$X, rep(0, nrow(d$X)), d$batch_vec, n_iter = 20, n_cores = n_cores)
  }
  a <- cont(2)
  b <- cont(2)
  expect_s3_class(a, "batchmix_fit_list")
  expect_length(a, 2)
  # 40 saved draws' worth (40 / 5 = 8) plus 20 / 5 = 4 new ones
  expect_equal(nrow(a[[1]]$samples), 12)
  expect_equal(vapply(a, function(ch) ch$Chain, numeric(1)), 1:2)
  for (i in 1:2) expect_identical(a[[i]]$samples, b[[i]]$samples)
  expect_false(identical(a[[1]]$samples, a[[2]]$samples))
  expect_false(is.null(attr(a, "convergence")))
})
