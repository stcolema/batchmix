#!/usr/bin/Rscript
#' @title Generate initial labels
#' @description For simulated data, generates an initial labelling for sampling.
#' @param alpha The mass in the stick breaking prior
#' @param K The number of classes available.
#' @param fixed The vector of 0s and 1s indicating which labels are to be held 
#' fixed.
#' @param labels The initial labelling. Defaults to NULL.
#' @return An N vector of labels.
#' @export
#' @examples
#' N <- 500
#' P <- 2
#' K <- 2
#' B <- 5
#' mean_dist <- 4
#' batch_dist <- 0.3
#' cluster_means <- seq(1, K) * mean_dist
#' batch_shift <- rnorm(B, mean = batch_dist, sd = batch_dist)
#' std_dev <- rep(2, K)
#' batch_var <- rep(1.2, B)
#' cluster_weights <- rep(1 / K, K)
#' batch_weights <- rep(1 / B, B)
#'
#' my_data <- generateBatchData(
#'   N,
#'   P,
#'   cluster_means,
#'   std_dev,
#'   batch_shift,
#'   batch_var,
#'   cluster_weights,
#'   batch_weights
#' )
#' 
#' initial_labels <- generateInitialLabels(1, K, my_data$fixed)
generateInitialLabels <- function(alpha, K, fixed, labels = NULL) {
  
  N <- length(fixed)
  if(is.null(labels)) {
    labels <- rep(0, N)
  }
  
  N_fixed <- sum(fixed)
  N_unfixed <- N - N_fixed
  labels_available <- seq(1, K)
  
  # Draw prior weights from the stick breaking prior
  w <- rStickBreakingPrior(alpha, K)
  
  labels[fixed == 0] <- sample(labels_available, N_unfixed, 
    replace = TRUE, 
    prob = w
  )
  
  labels
}
