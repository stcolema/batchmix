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
#' labels <- c(rep(1, 10), 
#'   sample(c(1,2), size = 40, replace = TRUE), 
#'   rep(2, 10), 
#'   sample(c(1,2), size = 40, replace = TRUE)
#' )
#' 
#' fixed <- c(rep(1, 10), rep(0, 40), rep(1, 10), rep(0, 40))
#' 
#' # Batch
#' batch_vec <- sample(seq(1, 5), replace = TRUE, size = 100)
#' 
#' # Sampling parameters
#' R <- 1000
#' burn <- 250
#' thin <- 50
#'
#' # MCMC samples
#' samples <- batchSemiSupervisedMixtureModel(X, R, thin, labels, fixed, batch_vec, "MVN")
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
  R <- mcmc_output$R
  thin <- mcmc_output$thin
  
  # Is the output semisupervised
  is_semisupervised <- mcmc_output$Semisupervised
  
  # What summary statistic is used to define our point estimates
  use_median <- point_estimate_method == "median"
  use_mean <- point_estimate_method == "mean"
  wrong_method <- ! (use_median | use_mean)
  if(wrong_method)
    stop("Wrong point estimate method given. Must be one of 'mean' or 'median'")
  
  # We burn the floor of burn / thin of these 
  eff_burn <- floor(burn / thin)
  
  # We record only the floor of R / thin samples
  eff_R <- floor(R / thin) - eff_burn
  
  # The indices dropped as part of the burn in
  dropped_indices <- seq(1, eff_burn)
  
  new_output <- mcmc_output
  
  # First apply a burn in to all quantities
  new_output$batch_corrected_data <- mcmc_output$batch_corrected_data[, , -dropped_indices]
  
  # The model fit measurements
  new_output$observed_likelihood <- mcmc_output$observed_likelihood[-dropped_indices, ]
  new_output$complete_likelihood <- mcmc_output$complete_likelihood[-dropped_indices, ]
  new_output$BIC <- mcmc_output$BIC[-dropped_indices, ]
  
  # The allocations and allocation probabilities
  new_output$samples <- mcmc_output$samples[-dropped_indices, ]
  
  if(is_semisupervised)
    new_output$alloc <- mcmc_output$alloc[, , -dropped_indices]
  
  # The sampled parameters
  new_output$means <- mcmc_output$means[, , -dropped_indices]
  new_output$covariance <- mcmc_output$covariance[ , , -dropped_indices]
  new_output$batch_shift <- mcmc_output$batch_shift[ , , -dropped_indices]
  new_output$batch_scale <- mcmc_output$batch_scale[ , , -dropped_indices]
  new_output$mean_sum <- mcmc_output$mean_sum[ , , -dropped_indices]
  new_output$cov_comb <- mcmc_output$cov_comb[ , , -dropped_indices]
  
  new_output$weights <- mcmc_output$weights[-dropped_indices, ]
  
  # The mean of the posterior samples for the parameters
  if(use_mean) {
    mean_est <- rowMeans(new_output$means, dims = 2L)
    shift_est <- rowMeans(new_output$batch_shift, dims = 2L)
    scale_est <- rowMeans(new_output$batch_scale, dims = 2L)
  }
  
  if(use_median) {
    mean_est <- apply(new_output$means, c(1, 2), stats::median)
    shift_est <- apply(new_output$batch_shift, c(1, 2), stats::median)
    scale_est <- apply(new_output$batch_scale, c(1, 2), stats::median)
  }
  
  if(type == "MVT") {
    new_output$t_df <- mcmc_output$t_df[-dropped_indices, ]
    if(use_mean) 
      new_output$t_df_est <- colMeans(new_output$t_df)
    if(use_median) 
      new_output$t_df_est <-apply(new_output$t_df, 2, stats::median)
  }
  
  # The covariance is represented as a matrix for reasons, but is more naturally 
  # thought of as a 3D array
  if(use_mean)
    cov_est <- rowMeans(new_output$covariance, dims = 2L)
  if(use_median)
    cov_est <- apply(new_output$covariance, c(1, 2), stats::median)
  
  cov_est_better_format <- array(0, c(P, P, K_max))
  
  # The indices for the columns corresponding to the first column for each 
  # clusters' covariance matrix, with one trailing index that is used as a bound
  cov_comb_cluster_inds <- cov_inds <- seq(1, P * (K_max + 1), by = P)
  
  for(k in cluster_inds) {
    lb <- cov_inds[k]
    ub <- cov_inds[k + 1] - 1
    columns_selected <- seq(lb, ub)
    cov_est_better_format[ , , k] <- cov_est[ , columns_selected]
  }
  
  # The combinations of the batch and group parameters are a little awkward
  # Do a nicer a format and record point estimates
  if(use_mean)
    mean_sum_est <- rowMeans(new_output$mean_sum, dims = 2L)
  if(use_median)
    mean_sum_est <- apply(new_output$mean_sum, c(1, 2), stats::median)
  
  mean_sum_better_format <- array(0, c(P, K_max, B))
  
  if(use_mean)
    cov_comb_est <- rowMeans(new_output$cov_comb, dims = 2L)
  if(use_median)
    cov_comb_est <- apply(new_output$cov_comb, c(1, 2), stats::median)
  
  cov_comb_better_format <- vector("list", B)
  cov_comb_better_format_entry <- array(0, c(P, P, K_max))
  
  cov_comb_batch_inds <- seq(0, (P + 1) * K_max * B, by = P * K_max)
  
  # Iterate over batches saving the mean sums and covariance combinations in a
  # more obvious fashion
  for(b in batch_inds) {

    mean_sum_better_format[, , b] <- mean_sum_est[, seq(b, b + K_max - 1)]
    
    rel_inds <- cov_comb_batch_inds[b] + cov_comb_cluster_inds
    for(k in cluster_inds) {
      lb <- rel_inds[k]
      ub <- rel_inds[k + 1] - 1
      columns_selected <- seq(lb, ub)
      
      cov_comb_better_format_entry[, , k] <- cov_comb_est[ , columns_selected]
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
  
  # The estimate of the inferred dataset
  if(use_mean)
    inferred_dataset <- rowMeans(new_output$batch_corrected_data, dims = 2L)
  if(use_median)
    inferred_dataset <- apply(new_output$batch_corrected_data, c(1, 2), stats::median)
  
  new_output$inferred_dataset <- inferred_dataset
  
  if(is_semisupervised) {
    # The estimate of the allocation probability matrix, the probability of the 
    # most probable class and the predicted class
    new_output$allocation_probability <- .alloc_prob <- calcAllocProb(new_output, 
      method = point_estimate_method
    )
    
    new_output$prob <- apply(.alloc_prob, 1, max)
    new_output$pred <- apply(.alloc_prob, 1, which.max)
  } 
  
  # Record the applied burn in
  new_output$burn <- burn 
  
  # Return the MCMC object with burn in applied and point estimates found
  new_output
}
