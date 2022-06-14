#!/usr/bin/Rscript
#' @title Generate initial labels
#' @description For simulated data, generates an initial labelling for sampling.
#' @param labels The true classes.
#' @param fixed The vector of 0s and 1s indicating which labels are to be held 
#' fixed.
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
#' initial_labels <- generateInitialLabels(my_data$group_IDs, my_data$fixed)
generateInitialLabels <- function(labels, fixed) {
  N <- length(labels)
  N_fixed <- sum(fixed)
  N_unfixed <- N - N_fixed
  ratio <- table(labels[fixed == 1]) / N_fixed
  labels[fixed == 0] <- sample(unique(labels), N_unfixed, 
    replace = TRUE, 
    prob = ratio
  )
  labels
}
