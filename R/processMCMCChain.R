#!/usr/bin/Rscript
#' @title Process MCMC chain
#' @description Applies a burn in to and finds a point estimate for the output
#' of ``batchSemiSupervisedMixtureModel``.
#' @param mcmc_output Output from ``batchSemiSupervisedMixtureModel``
#' @param burn The number of MCMC samples to drop as part of a burn in.
#' @param point_estimate_method Summary statistic used to define the point
#' estimate. Must be ``'mean'`` or ``'median'``. ``'median'`` is the default.
#' @returns A named list similar to the output of
#' ``batchSemiSupervisedMixtureModel`` with some additional entries:
#'
#'  * ``mean_est``: $(P x K)$ matrix. The point estimate of the cluster
#'  means with columns  corresponding to clusters.
#'
#'  * ``cov_est``: $(P x P x K)$ array. The point estimate of the
#'  cluster covariance matrices with slices corresponding to clusters.
#'
#'  * ``shift_est``: $(P x B)$ matrix. The point estimate of the batch
#'  shift effect with columns  corresponding to batches.
#'
#'  * ``scale_est``: $(P x B)$ matrix. The point estimate of the batch
#'  scale effects. The $bth$ column contains the diagonal entries of the scaling
#'  matrix for the $bth£ batch.
#'
#'  * ``mean_sum_est``: $(P x K x B)$ array. The point estimate of the
#'  sum of the cluster  means and the batch shift effect with columns
#'  corresponding to clusters and slices to batches.
#'
#'  * ``cov_comb_est``: List of length $B$, with each entry being a
#'  $(P x P x K)$ array. The point estimate of the combination of the
#'  cluster covariance matrices and the batch scale effect with list entries
#'  corresponding to batches and slices of each array corresponding to clusters.
#'
#'  * ``inferred_dataset``: $(N x P)$ matrix. The inferred ``batch-free''
#'  dataset.
#'
#'  * ``allocation_probability``: $(N x K)$ matrix. The point estimate of
#'  the allocation probabilities for each data point to each class.
#'
#'  * ``prob``: $N$ vector. The point estimate of the probability of being
#'  allocated to the class with the highest probability.
#'
#'  * ``pred``: $N$ vector. The predicted class for each sample.
#'
#' @export
#' @examples
#'
#' # Data in a matrix format
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
#' burn <- 250
#' thin <- 50
#'
#' # MCMC samples
#' samples <- runBatchMix(X, n_iter, thin, batch_vec, "MVN",
#'   initial_labels = labels,
#'   fixed = fixed
#' )
#'
#' # Process the MCMC samples
#' processed_samples <- processMCMCChain(samples, burn)
#'
#' @importFrom stats median
processMCMCChain <- function(mcmc_output, burn, point_estimate_method = "median") {
  # Dimensions of the dataset
  N <- mcmc_output$N
  P <- mcmc_output$P
  K_max <- mcmc_output$K_max
  B <- mcmc_output$B

  # The type of mixture model used
  type <- mcmc_output$type

  # Indices for clusters and batches
  batch_inds <- seq(1, B)
  cluster_inds <- seq(1, K_max)

  # MCMC iterations and thinning
  n_iter <- mcmc_output$n_iter
  thin <- mcmc_output$thin

  # What summary statistic is used to define our point estimates
  use_median <- point_estimate_method == "median"
  use_mean <- point_estimate_method == "mean"
  wrong_method <- !(use_median | use_mean)
  if (wrong_method) {
    stop("Wrong point estimate method given. Must be one of 'mean' or 'median'")
  }

  # We burn the floor of burn / thin of these
  eff_burn <- floor(burn / thin)

  # We record only the floor of n_iter / thin samples
  eff_R <- floor(n_iter / thin) - eff_burn

  # The indices to KEEP after burn-in. This is deliberately not
  # `-seq_len(eff_burn)` used directly at each call site: when eff_burn is 0
  # (e.g. any 0 < burn < thin), `-seq_len(0)` is `-integer(0)`, and indexing
  # with an empty vector - negated or not - selects NOTHING in R, not
  # "everything" (`(1:5)[-integer(0)]` is `integer(0)`, not `1:5`) - so using
  # it directly would silently drop every sample instead of none. TRUE
  # (recycled, i.e. "keep everything") sidesteps this for eff_burn == 0;
  # -seq_len(eff_burn) is used otherwise.
  keep_indices <- if (eff_burn > 0) -seq_len(eff_burn) else TRUE

  new_output <- mcmc_output

  # First apply a burn in to all quantities
  # drop = FALSE matters here: without it, R silently collapses the result
  # to a 2D matrix whenever P == 1 (a perfectly ordinary univariate dataset)
  # since indexing drops any resulting dimension of extent 1, not just the
  # one actually being subset - this broke rowMeans(..., dims = 2L) below
  # outright for point_estimate_method = "mean", and silently returned the
  # wrong shape (no averaging across iterations at all) for the median path.
  new_output$batch_corrected_data <- mcmc_output$batch_corrected_data[, , keep_indices, drop = FALSE]

  # The model fit measurements
  new_output$observed_likelihood <- mcmc_output$observed_likelihood[keep_indices, ]
  new_output$complete_likelihood <- mcmc_output$complete_likelihood[keep_indices, ]
  new_output$BIC <- mcmc_output$BIC[keep_indices, ]
  
  # The allocations and allocation probabilities. alloc is populated for
  # every unfixed item regardless of Semisupervised (a fully unsupervised
  # fit has every item unfixed), so it is trimmed unconditionally - not
  # doing so left it at its full, un-burned length for unsupervised fits.
  new_output$samples <- mcmc_output$samples[keep_indices, ]
  # drop = FALSE for the same reason as batch_corrected_data above (this
  # collapses whenever K == 1).
  new_output$alloc <- mcmc_output$alloc[, , keep_indices, drop = FALSE]

  # The sampled parameters
  new_output$means <- mcmc_output$means[, , keep_indices, drop = FALSE]
  new_output$covariance <- mcmc_output$covariance[, , keep_indices, drop = FALSE]
  new_output$batch_shift <- mcmc_output$batch_shift[, , keep_indices, drop = FALSE]
  new_output$batch_scale <- mcmc_output$batch_scale[, , keep_indices, drop = FALSE]
  new_output$mean_sum <- mcmc_output$mean_sum[, , keep_indices, drop = FALSE]
  new_output$cov_comb <- mcmc_output$cov_comb[, , keep_indices, drop = FALSE]

  new_output$weights <- mcmc_output$weights[keep_indices, , drop = FALSE]

  m_scale_sampled <- mcmc_output$sample_m_scale
  if(m_scale_sampled) {
    new_output$lambda_2 <- mcmc_output$lambda_2[keep_indices]
  }

  # Batch x cluster interaction term and GP-correlated batch weights (if
  # requested when the chain was run) need the same burn-in applied as
  # every other sampled quantity above.
  interaction_used <- isTRUE(mcmc_output$include_interaction)
  if (interaction_used) {
    new_output$gamma <- mcmc_output$gamma[, , keep_indices, drop = FALSE]
  }

  correlated_weights_used <- isTRUE(mcmc_output$weight_prior_type > 0)
  if (correlated_weights_used) {
    new_output$w_batch <- mcmc_output$w_batch[, , keep_indices, drop = FALSE]
    new_output$eta_alr <- mcmc_output$eta_alr[, , keep_indices, drop = FALSE]
    new_output$gp_tau2 <- mcmc_output$gp_tau2[keep_indices]
    new_output$gp_length_scale <- mcmc_output$gp_length_scale[keep_indices]
    new_output$pp_mu <- mcmc_output$pp_mu[, keep_indices, drop = FALSE]
    new_output$pp_tau2 <- mcmc_output$pp_tau2[, keep_indices, drop = FALSE]
    new_output$gp_beta <- mcmc_output$gp_beta[, keep_indices, drop = FALSE]
  }

  if (type == "MVT") {
    # Trimmed (but not yet relabelled - see below) here, alongside every
    # other cluster-indexed quantity, so relabelChain() picks it up too;
    # trimming it after relabelChain() would silently leave t_df/t_df_est
    # exposed to the same label-switching bug this whole block exists to fix.
    new_output$t_df <- mcmc_output$t_df[keep_indices, , drop = FALSE]
  }

  # Mixture models are only identified up to a permutation of the component
  # labels ("label switching"). Averaging/taking the median of the raw,
  # cluster-indexed arrays below is only valid once every iteration's labels
  # have been aligned to a common reference - see relabelChain() for the
  # method and why this must happen before, not after, the point estimates
  # computed from here on.
  new_output <- relabelChain(new_output, K_max = K_max)

  # The mean of the posterior samples for the parameters
  if (use_mean) {
    mean_est <- rowMeans(new_output$means, dims = 2L)
    shift_est <- rowMeans(new_output$batch_shift, dims = 2L)
    scale_est <- rowMeans(new_output$batch_scale, dims = 2L)
  }

  if (use_median) {
    mean_est <- apply(new_output$means, c(1, 2), stats::median)
    shift_est <- apply(new_output$batch_shift, c(1, 2), stats::median)
    scale_est <- apply(new_output$batch_scale, c(1, 2), stats::median)
  }

  if (type == "MVT") {
    if (use_mean) {
      new_output$t_df_est <- colMeans(new_output$t_df)
    }
    if (use_median) {
      new_output$t_df_est <- apply(new_output$t_df, 2, stats::median)
    }
  }

  # The covariance is represented as a matrix for reasons, but is more naturally
  # thought of as a 3D array
  if (use_mean) {
    cov_est <- rowMeans(new_output$covariance, dims = 2L)
  }
  if (use_median) {
    cov_est <- apply(new_output$covariance, c(1, 2), stats::median)
  }

  cov_est_better_format <- array(0, c(P, P, K_max))

  # The indices for the columns corresponding to the first column for each
  # clusters' covariance matrix, with one trailing index that is used as a bound
  cov_inds <- seq(1, P * (K_max + 1), by = P)

  for (k in cluster_inds) {
    lb <- cov_inds[k]
    ub <- cov_inds[k + 1] - 1
    columns_selected <- seq(lb, ub)
    cov_est_better_format[, , k] <- cov_est[, columns_selected]
  }

  # The combinations of the batch and group parameters are a little awkward
  # Do a nicer a format and record point estimates
  if (use_mean) {
    mean_sum_est <- rowMeans(new_output$mean_sum, dims = 2L)
  }
  if (use_median) {
    mean_sum_est <- apply(new_output$mean_sum, c(1, 2), stats::median)
  }

  mean_sum_better_format <- array(0, c(P, K_max, B))

  if (use_mean) {
    cov_comb_est <- rowMeans(new_output$cov_comb, dims = 2L)
  }
  if (use_median) {
    cov_comb_est <- apply(new_output$cov_comb, c(1, 2), stats::median)
  }

  cov_comb_better_format <- vector("list", B)

  # mean_sum/cov_comb are stored with the cluster index varying SLOWEST (the
  # C++ sampler indexes both by kb = k * B + b - see matrixCombinations() in
  # mvnSampler.cpp/mvnSamplerSeparationStrategy.cpp), i.e. for a fixed batch
  # b, the K_max cluster blocks are strided by B blocks apart, NOT
  # contiguous. mean_sum's blocks are 1 column wide; cov_comb's are P
  # columns wide.
  for (b in batch_inds) {
    mean_sum_cols <- seq(b, by = B, length.out = K_max)
    mean_sum_better_format[, , b] <- mean_sum_est[, mean_sum_cols, drop = FALSE]

    cov_comb_better_format_entry <- array(0, c(P, P, K_max))
    for (k in cluster_inds) {
      block <- (k - 1) * B + (b - 1) # 0-indexed block number, matching k * B + b in the sampler
      lb <- block * P + 1
      ub <- lb + P - 1
      cov_comb_better_format_entry[, , k] <- cov_comb_est[, lb:ub]
    }
    cov_comb_better_format[[b]] <- cov_comb_better_format_entry
  }

  # Save hte estimated parameters to the output object
  new_output$mean_est <- mean_est
  new_output$shift_est <- shift_est
  new_output$scale_est <- scale_est
  new_output$cov_est <- cov_est_better_format
  new_output$mean_sum_est <- mean_sum_better_format
  new_output$cov_comb_est <- cov_comb_better_format

  if (interaction_used) {
    if (use_mean) {
      new_output$gamma_est <- rowMeans(new_output$gamma, dims = 2L)
    } else {
      new_output$gamma_est <- apply(new_output$gamma, c(1, 2), stats::median)
    }
  }

  if (correlated_weights_used) {
    if (use_mean) {
      new_output$w_batch_est <- rowMeans(new_output$w_batch, dims = 2L)
    } else {
      new_output$w_batch_est <- apply(new_output$w_batch, c(1, 2), stats::median)
    }
  }

  # The estimate of the inferred dataset
  if (use_mean) {
    inferred_dataset <- rowMeans(new_output$batch_corrected_data, dims = 2L)
  }
  if (use_median) {
    inferred_dataset <- apply(new_output$batch_corrected_data, c(1, 2), stats::median)
  }

  new_output$inferred_dataset <- inferred_dataset

  # The estimate of the allocation probability matrix, the probability of
  # the most probable class and the predicted class. alloc is populated for
  # every unfixed item regardless of Semisupervised (see above), and is
  # already relabelled (relabelChain() ran on new_output before any point
  # estimate here), so this is just as well-defined - and just as useful,
  # e.g. for comparing predicted labels against a known ground truth on
  # simulated data - for a fully unsupervised fit as a semi-supervised one.
  new_output$allocation_probability <- .alloc_prob <- calcAllocProb(new_output,
    method = point_estimate_method
  )

  new_output$prob <- apply(.alloc_prob, 1, max)
  new_output$pred <- apply(.alloc_prob, 1, which.max)

  # Record the applied burn in
  new_output$burn <- burn

  # Return the MCMC object with burn in applied and point estimates found
  new_output
}
