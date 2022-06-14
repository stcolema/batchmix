#!/usr/bin/Rscript
#
#' @title Generate group IDs
#' @description Generate group IDs within ``generateBatchData``.
#' @param N The number of items (rows) to generate.
#' @param K The number of groups to genetare.
#' @param B The number of batches present in ``batch_IDs``.
#' @param batch_IDs The batch membership of each item.
#' @param group_weights One of either a K x B matrix of the expected proportion 
#' of each batch in each group or a K-vector of the expected proportion of the
#' entire dataset in each group.
#' @param varying_group_within_batch Flag indicating if the groups are vvarying 
#' across batches.
#' @return A N-vector of group membership.
#' @examples
#' N <- 500
#' K <- 2
#' B <- 5
#' group_weights <- rep(1 / K, K)
#' batch_weights <- rep(1 / B, B)
#' batch_IDs <- sample(seq(1, B), N, replace = TRUE, prob = batch_weights)
#' varying_group_within_batch <- FALSE
#' group_IDs <- generateGroupIDsInSimulator(
#'   N,
#'   K,
#'   B,
#'   batch_IDs, 
#'   group_weights, 
#'   varying_group_within_batch
#' )
generateGroupIDsInSimulator <- function(N,
                                        K, 
                                        B, 
                                        batch_IDs, 
                                        group_weights, 
                                        varying_group_within_batch) {
  # Generate group membership, potentially allowing imbalance across batches
  group_IDs <- rep(0, N)
  if (varying_group_within_batch) {
    # Iterate over the batches to sample appropriate labels
    for (b in seq(1, B)) {
      batch_ind <- which(batch_IDs == b)
      N_b <- length(batch_ind)
      
      group_IDs[batch_ind] <- sample(seq(1, K), N_b,
        replace = TRUE,
        prob = group_weights[, b]
      )
    }
  } else {
    # The membership vector for the N points
    group_IDs <- sample(seq(1, K), N, replace = TRUE, prob = group_weights)
  }
  group_IDs
}