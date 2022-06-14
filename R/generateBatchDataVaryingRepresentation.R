#!/usr/bin/Rscript
#' @title Generate batch data
#' @description Generate data from groups across batches. Assumes independence
#' across columns. In each column the parameters are randomly permuted for both
#' the groups and batches.
#' @param N The number of items (rows) to generate.
#' @param P The number of columns in the generated dataset.
#' @param group_means A vector of the group means for a column.
#' @param group_std_dev A vector of group standard deviations for a column.
#' @param batch_shift A vector of batch means in a column.
#' @param batch_scale A vector of batch standard deviations within a column.
#' @param group_weights A K x B matrix of the expected proportion of N in each group in each batch.
#' @param batch_weights A vector of the expected proportion of N in each batch.
#' @param frac_known The expected fraction of observed labels. Used to generate
#' a ``fixed`` vector to feed into the ``batchSemiSupervisedMixtureModel`` function.
#' @return A list of 4 objects; the data generated from the groups with and
#' without batch effects, the label indicating the generating group and the
#' batch label.
#' @export
#' @examples
#' N <- 500
#' P <- 2
#' K <- 2
#' B <- 5
#' mean_dist <- 4
#' batch_dist <- 0.3
#' group_means <- seq(1, K) * mean_dist
#' batch_shift <- rnorm(B, mean = batch_dist, sd = batch_dist)
#' std_dev <- rep(2, K)
#' batch_var <- rep(1.2, B)
#' group_weights <- matrix(c(
#'   0.8, 0.6, 0.4, 0.2, 0.2,
#'   0.2, 0.4, 0.6, 0.8, 0.8
#' ),
#' nrow = K, ncol = B, byrow = TRUE
#' )
#' batch_weights <- rep(1 / B, B)
#'
#' my_data <- generateBatchDataVaryingRepresentation(
#'   N,
#'   P,
#'   group_means,
#'   std_dev,
#'   batch_shift,
#'   batch_var,
#'   group_weights,
#'   batch_weights
#' )
#' @importFrom stats rnorm
generateBatchDataVaryingRepresentation <- function(N,
                                                   P,
                                                   group_means,
                                                   group_std_dev,
                                                   batch_shift,
                                                   batch_scale,
                                                   group_weights,
                                                   batch_weights,
                                                   frac_known = 0.2) {

  # The number of groups to generate
  K <- length(group_means)

  # The number of batches to generate
  B <- length(batch_shift)

  if (ncol(group_weights) != B) {
    stop("Number of columns in group weight matrix does not match the number of batches.")
  }

  if (nrow(group_weights) != K) {
    stop("Number of rows in group weight matrix does not match the number of groupes.")
  }

  # The membership vector for the N points, currently empty
  labels <- rep(0, N)

  # The batch labels for the N points
  batches <- sample(seq(1, B), N, replace = TRUE, prob = batch_weights)

  # The fixed labels for the semi-supervised case
  fixed <- sample(seq(0, 1), N, 
    replace = TRUE, 
    prob = c(1 - frac_known, frac_known)
  )

  # The data matrices
  true_data <- observed_data <- matrix(0, nrow = N, ncol = P)

  # Iterate over the batches to sample appropriate labels
  for (b in seq(1, B)) {
    batch_ind <- which(batches == b)
    N_b <- length(batch_ind)
    
    labels[batch_ind] <- sample(seq(1, K), N_b, 
      replace = TRUE, 
      prob = group_weights[, b]
    )
    
  }

  # Generate the data
  for (p in seq(1, P)) {

    # To provide different information in each column, randomly sample the
    # parameters with each group and batch
    reordered_group_means <- sample(group_means)
    reordered_group_std_devs <- sample(group_std_dev)

    reordered_batch_shift <- sample(batch_shift)
    reordered_batch_scale <- sample(batch_scale)

    for (n in seq(1, N)) {

      # Find the group and batch
      b <- batches[n]
      k <- labels[n]

      # Random data
      x <- stats::rnorm(1)

      # For ease of reading the following lines, create group and batch parameters
      .mu <- reordered_group_means[k]
      .sd <- reordered_group_std_devs[k]
      .m <- reordered_batch_shift[b]
      .s <- reordered_batch_scale[b]
      
      # Corrected and observed data point
      true_data[n, p] <- x * .sd + .mu
      observed_data[n, p] <- x * .sd * .s + .mu + .m
    }
  }

  # Return a list of the group labels, batch labels, fixed points and the
  # observed and true datasets.
  list(
    observed_data = observed_data,
    corrected_data = true_data,
    group_IDs = labels,
    batch_IDs = batches,
    fixed = fixed
  )
}
