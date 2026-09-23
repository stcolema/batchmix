#!/usr/bin/Rscript
#' @title Calculate allocation probabilities
#' @description Calculate the empirical allocation probability for each class
#' based on the sampled allocation probabilities.
#' @param mcmc_samples Output from ``batchSemiSupervisedMixtureModel``.
#' @param burn The number of samples to discard.
#' @param method The point estimate to use. ``method = 'mean'`` or
#' ``method = 'median'``. ``'median'`` is the default.
#' @return An N x K matrix of class probabilities.
#' @export
#' @examples
#'
#' # Data in matrix format
#' X <- matrix(c(rnorm(100, 0, 1), rnorm(100, 3, 1)), ncol = 2, byrow = TRUE)
#'
#' # Initial labelling
#' labels <- c(
#'   rep(1, 10),
#'   sample(c(1, 2), size = 40, replace = TRUE),
#'   rep(2, 10),
#'   sample(c(1, 2), size = 40, replace = TRUE)
#' )
#'
#' fixed <- c(rep(1, 10), rep(0, 40), rep(1, 10), rep(0, 40))
#'
#' # Batch
#' batch_vec <- sample(seq(1, 5), replace = TRUE, size = 100)
#'
#' # Sampling parameters
#' n_iter <- 1000
#' thin <- 50
#'
#' # MCMC samples and BIC vector
#' samples <- batchSemiSupervisedMixtureModel(X, n_iter, thin, labels, fixed, batch_vec, "MVN")
#'
#' # Burn in
#' burn <- 20
#' eff_burn <- burn / thin
#'
#' # Probability across classes
#' probs <- calcAllocProb(samples, burn = burn)
#'
calcAllocProb <- function(mcmc_samples, burn = 0, method = "median") {
  n_iter <- mcmc_samples$n_iter
  thin <- mcmc_samples$thin
  .alloc <- mcmc_samples$alloc
  .samples <- mcmc_samples$samples

  if (burn > 0) {
    if (burn > n_iter) {
      stop("Burn in exceeds number of iterations run.")
    }

    eff_burn <- floor(burn / thin)
    # Deliberately not `-seq_len(eff_burn)` used directly: when eff_burn is
    # 0 (e.g. 0 < burn < thin), indexing with an empty vector - negated or
    # not - selects NOTHING in R, not "everything" (`(1:5)[-integer(0)]` is
    # `integer(0)`, not `1:5`), so that would silently drop every sample
    # instead of none. TRUE (recycled) keeps everything for eff_burn == 0.
    keep_samples <- if (eff_burn > 0) -seq_len(eff_burn) else TRUE
    # drop = FALSE: without it, indexing collapses .alloc to a 2D matrix
    # whenever K == 1, breaking every array-shaped operation below.
    .alloc <- .alloc[, , keep_samples, drop = FALSE]
    .samples <- .samples[keep_samples, , drop = FALSE]
  }

  # Mixture models are only identified up to a permutation of the component
  # labels ("label switching") - align every iteration's allocation
  # probabilities to a common reference before averaging/taking their
  # median, or the summary below conflates genuine allocation uncertainty
  # with harmless label permutation. See relabelChain(); calling it here
  # with only `samples`/`alloc` set is idempotent (a no-op) when `.alloc`
  # has already been relabelled by a caller such as processMCMCChain().
  K_max <- if (!is.null(mcmc_samples$K_max)) mcmc_samples$K_max else ncol(.alloc)
  .alloc <- relabelChain(list(samples = .samples, alloc = .alloc, K_max = K_max))$alloc

  probs <- NULL

  if (method == "median") {
    probs <- apply(.alloc, c(1, 2), median)
  }
  if (method == "mean") {
    probs <- rowSums(.alloc, dims = 2) / dim(.alloc)[3]
  }
  if (length(probs) == 1) {
    if (is.null(probs)) {
      stop("``method`` must be one of 'mean' or 'median'")
    }
  }
  probs
}
