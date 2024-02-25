#' @title Predict from multiple MCMC chains
#' @description Applies a burn in to and finds a point estimate by combining
#' multiple chains of ``callMDI``.
#' @param mcmc_outputs Output from ``runMCMCChains``
#' @param burn The number of MCMC samples to drop as part of a burn in.
#' @param point_estimate_method Summary statistic used to define the point
#' estimate. Must be ``'mean'`` or ``'median'``. ``'median'`` is the default.
#' @param chains_already_processed Logical indicating if the chains have already
#' had a burn-in applied.
#' @returns A named list of quantities related to prediction/clustering:
#'
#'  * ``allocation_probability``: List with an $(N x K)$ matrix if the model is
#'    semi-supervised. The point estimate of the allocation probabilities for
#'    each data point to each class.
#'
#'  * ``prob``: $N$ vector of the point estimate of the probability of being
#'    allocated to the class with the highest probability.
#'
#'  * ``pred``: $N$ vector of the predicted class for each sample. If the model
#'    is unsupervised then the ``salso`` function from Dahl et al. (2021) is
#'    used on the sampled partitions using the default settings.
#'
#'  * ``samples``: List of sampled allocations for each view. Columns
#'    correspond to items being clustered, rows to MCMC samples.
#' @examples
#' \donttest{
#' # Data dimensions
#' N <- 600
#' P <- 4
#' K <- 5
#' B <- 7
#'
#' # Generating model parameters
#' mean_dist <- 2.25
#' batch_dist <- 0.3
#' group_means <- seq(1, K) * mean_dist
#' batch_shift <- rnorm(B, mean = batch_dist, sd = batch_dist)
#' std_dev <- rep(2, K)
#' batch_var <- rep(1.2, B)
#' group_weights <- rep(1 / K, K)
#' batch_weights <- rep(1 / B, B)
#' dfs <- c(4, 7, 15, 60, 120)
#'
#' my_data <- generateBatchData(
#'   N,
#'   P,
#'   group_means,
#'   std_dev,
#'   batch_shift,
#'   batch_var,
#'   group_weights,
#'   batch_weights,
#'   type = "MVT",
#'   group_dfs = dfs
#' )
#'
#'
#' X <- my_data$observed_data
#'
#' true_labels <- my_data$group_IDs
#' fixed <- my_data$fixed
#' batch_vec <- my_data$batch_IDs
#'
#' alpha <- 1
#' initial_labels <- generateInitialLabels(alpha, K, fixed, true_labels)
#'
#' # Sampling parameters
#' R <- 1000
#' thin <- 25
#' burn <- 100
#' n_chains <- 2
#'
#' # Density choice
#' type <- "MVT"
#'
#' # MCMC samples and BIC vector
#' mcmc_outputs <- runMCMCChains(
#'   X,
#'   n_chains,
#'   R,
#'   thin,
#'   batch_vec,
#'   type,
#'   initial_labels = initial_labels,
#'   fixed = fixed
#' )
#' ensemble_mod <- predictFromMultipleChains(mcmc_outputs, burn)
#' }
#' @export
predictFromMultipleChains <- function(mcmc_outputs,
                                      burn,
                                      point_estimate_method = "median",
                                      chains_already_processed = FALSE) {
  if (chains_already_processed) {
    processed_chains <- mcmc_outputs
  } else {
    # Process the chains, making point estimates and applying a burn-in
    processed_chains <- processMCMCChains(mcmc_outputs, burn, point_estimate_method)
  }

  # The number of chains
  n_chains <- length(processed_chains)
  chain_indices <- seq(1, n_chains)

  # This is used to derive some characteristics of the model run
  first_chain <- processed_chains[[1]]

  # Dimensions of the dataset
  N <- first_chain$N
  K <- first_chain$K

  # MCMC call
  R <- first_chain$R
  thin <- first_chain$thin

  # The type of mixture model used
  type <- first_chain$type

  # Is the output semisupervised
  is_semisupervised <- first_chain$Semisupervised

  # What summary statistic is used to define our point estimates
  use_median <- point_estimate_method == "median"
  use_mean <- point_estimate_method == "mean"
  wrong_method <- !(use_median | use_mean)
  if (wrong_method) {
    stop("Wrong point estimate method given. Must be one of 'mean' or 'median'")
  }

  # We burn the floor of burn / thin of these
  eff_burn <- floor(burn / thin)

  # We record only the floor of R / thin samples
  eff_R <- floor(R / thin) - eff_burn

  # The indices dropped as part of the burn in
  dropped_indices <- seq(1, eff_burn)

  # Setup the output list
  merged_outputs <- list()

  merged_outputs$R <- R
  merged_outputs$thin <- thin
  merged_outputs$burn <- burn
  merged_outputs$n_chains <- n_chains
  merged_outputs$Semisupervised <- is_semisupervised

  merged_outputs$Point_estimate <- point_estimate_method

  merged_outputs$N <- N
  merged_outputs$K <- K

  merged_outputs$weights <- do.call(rbind, lapply(processed_chains, function(x) x$weights))
  merged_outputs$concentration <- do.call(rbind, lapply(processed_chains, function(x) x$concentration))

  first_chain <- TRUE

  merged_outputs$allocation_probability <- .alloc_prob <- matrix(
    0,
    N,
    K
  )
  for (ii in chain_indices) {
    .curr_chain <- processed_chains[[ii]]
    in_first_chain <- (ii == 1)

    if (in_first_chain) {
      merged_outputs$samples <- .curr_chain$samples[, , drop = TRUE]
    } else {
      .prev <- merged_outputs$samples
      .current <- .curr_chain$samples[, , drop = TRUE]
      merged_outputs$samples <- rbind(.prev, .current)
    }

    if (is_semisupervised) {
      .prev <- .alloc_prob
      .curr <- .curr_chain$allocation_probability

      .alloc_prob <- .prev + .curr
    }
  }

  if (is_semisupervised) {
    # Normalise the probabilities
    .alloc_prob <- .alloc_prob / n_chains

    merged_outputs$allocation_probability <- .alloc_prob

    merged_outputs$prob <- .prob <- apply(.alloc_prob, 1, max)
    merged_outputs$pred <- apply(.alloc_prob, 1, which.max)
  } else {
    merged_outputs$psm <- .psm <- createSimilarityMat(merged_outputs$samples)
    merged_outputs$pred <- minVI(.psm, merged_outputs$samples, method = "avg") # suppressWarnings(salso::salso(merged_outputs$samples))
  }
  merged_outputs
}
