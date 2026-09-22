# Deterministic, small inputs shared by the characterization tests. Kept
# separate from the package's own data generators (generateBatchData() etc.)
# so that these fixtures stay stable even if the generators change - the
# characterization tests should only be sensitive to sampler behaviour.

make_characterization_data <- function(seed = 20240101, N = 45, P = 2, K = 3, B = 2) {
  set.seed(seed)

  X <- matrix(rnorm(N * P, mean = rep(c(0, 3), length.out = P)), nrow = N, ncol = P)
  batch_vec <- sample(seq_len(B), N, replace = TRUE)
  labels <- sample(seq_len(K), N, replace = TRUE)

  # A reproducible ~30% of items fixed for semi-supervised scenarios.
  fixed_semi <- as.integer(seq_len(N) %% 3 == 0)
  fixed_none <- rep(0L, N)

  list(
    N = N, P = P, K = K, B = B,
    X = X,
    batch_vec = batch_vec,
    labels = labels,
    fixed_semi = fixed_semi,
    fixed_none = fixed_none
  )
}

# A small mixed-type dataset (one continuous column, one binary column, with
# a couple of missing entries) for the MVN_MIXED family.
make_mixed_characterization_data <- function(seed = 20240101, N = 45, K = 3, B = 2) {
  d <- make_characterization_data(seed = seed, N = N, P = 2, K = K, B = B)

  X <- d$X
  X[, 2] <- as.numeric(X[, 2] > median(X[, 2]))

  # A handful of missing entries, deterministic given the seed.
  set.seed(seed + 1)
  na_idx <- sample(seq_len(N), size = 3)
  X[na_idx, 1] <- NA_real_

  d$X <- X
  d$column_type <- c(0L, 1L)
  d$censor_code <- matrix(0L, nrow = N, ncol = 2)
  d
}
