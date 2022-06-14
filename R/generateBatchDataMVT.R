#!/usr/bin/Rscript
#
#' @title Generate batch data from a multivariate t distribution
#' @description Generate data from K multivariate t distributions with 
#' additional noise from batches. Assumes independence across columns. In each
#' column the parameters are randomly permuted for both the groups and batches.
#' @param N The number of items (rows) to generate.
#' @param P The number of columns in the generated dataset.
#' @param group_means A vector of the group means for a column.
#' @param group_std_devs A vector of group standard deviations for a column.
#' @param batch_shift A vector of batch means in a column.
#' @param batch_scale A vector of batch standard deviations within a column.
#' @param group_weights A K x B matrix of the expected proportion of N in each group in each batch.
#' @param batch_weights A vector of the expected proportion of N in each batch.
#' @param frac_known The number of items with known labels.
#' @param dfs A K-vector of the group specific degrees of freedom.
#' @return A list of 5 objects; the data generated from the groups with and
#' without batch effects, the label indicating the generating group, the
#' batch label and the vector indicating training versus test.
#' @export
#' @examples
#' N <- 500
#' P <- 2
#' K <- 2
#' B <- 5
#' mean_dist <- 4
#' batch_dist <- 0.3
#' group_means <- seq(1, K)* mean_dist
#' batch_shift <- rnorm(B, mean = batch_dist, sd = batch_dist)
#' std_dev <- rep(2, K)
#' batch_var <- rep(1.2, B)
#' group_weights <- rep(1 / K, K)
#' batch_weights <- rep(1 / B, B)
#' dfs <- c(4, 7)
#' my_data <- generateBatchDataMVT(
#'   N,
#'   P,
#'   group_means,
#'   std_dev,
#'   batch_shift,
#'   batch_var,
#'   group_weights,
#'   batch_weights,
#'   dfs
#' )
#' @importFrom stats rnorm rchisq
generateBatchDataMVT <- function(N,
                                 P,
                                 group_means,
                                 group_std_devs,
                                 batch_shift,
                                 batch_scale,
                                 group_weights,
                                 batch_weights,
                                 dfs,
                                 frac_known = 0.2) {
  
  # The number of groups to generate
  K <- length(group_means)
  
  # The number of batches to generate
  B <- length(batch_shift)
  
  # The membership vector for the N points
  group_IDs <- sample(seq(1, K), N, replace = TRUE, prob = group_weights)
  
  # The batch labels for the N points
  batch_IDs <- sample(seq(1, B), N, replace = TRUE, prob = batch_weights)
  
  # Fixed labels
  fixed <- sample(seq(0, 1), N, 
    replace = TRUE, 
    prob = c(1 - frac_known, frac_known)
  )
  
  # The data matrices
  observed_data <- true_data <- matrix(nrow = N, ncol = P)
  
  # Iterate over each of the columns permuting the means associated with each
  # label.
  for (p in seq(1, P))
  {
    
    # To provide different information in each column, randomly sample the 
    # parameters with each group and batch
    reordered_group_means <- sample(group_means)
    reordered_group_std_devs <- sample(group_std_devs)
    
    reordered_batch_shift <- sample(batch_shift)
    reordered_batch_scale <- sample(batch_scale)
    
    # Draw n points from the K univariate Gaussians defined by the permuted means.
    for (n in seq(1, N)) {
      
      # Draw a point from a standard normal
      x <- stats::rnorm(1)
      k <- group_IDs[n]
      b <- batch_IDs[n]
      
      chi_draw <- stats::rchisq(1, dfs[k])
      
      # For ease of reading the following lines, create group and batch parameters
      .mu <- reordered_group_means[k]
      .sd <- reordered_group_std_devs[k]
      .m <- reordered_batch_shift[b]
      .s <- reordered_batch_scale[b]
      
      # Adjust to the group distribution
      true_data[n, p] <- x * .sd * sqrt(dfs[k] / chi_draw) + .mu
      
      # Adjust to the batched group distribution
      observed_data[n, p] <- x * .sd * .s * sqrt(dfs[k] / chi_draw) + .mu + .m
    }
  }

  # Return the data, the data without batch effects, the allocation labels and
  # the batch labels.
  list(
    observed_data = observed_data,
    corrected_data = true_data,
    group_IDs = group_IDs,
    batch_IDs = batch_IDs,
    fixed = fixed
  )
}

