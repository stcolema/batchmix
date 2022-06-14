#!/usr/bin/Rscript
#
#' @title Check data generation inputs
#' @description Checks that the inputs for the ``generateBatchData`` function 
#' are correct. For internal use only.
#' @param N The number of items (rows) to generate.
#' @param P The number of columns in the generated dataset.
#' @param group_means A vector of the group means for a column.
#' @param group_std_devs A vector of group standard deviations for a column.
#' @param batch_shift A vector of batch means in a column.
#' @param batch_scale A vector of batch standard deviations within a column.
#' @param group_weights A K x B matrix of the expected proportion of N in each group in each batch.
#' @param batch_weights A vector of the expected proportion of N in each batch.
#' @param type A string indicating if data should be generated from
#' multivariate normal ("MVN") or multivariate t ("MVT") densities.
#' @param group_dfs A K-vector of the group specific degrees of freedom.
#' @param frac_known The number of items with known labels.
#' @param permute_variables Logical indicating if group and batch means and
#' standard deviations should be permuted in each column or not.
#' @param scale_data Logical indicating if data should be mean centred and 
#' standardised.
#' @return NULL
#' @examples
#' N <- 500
#' P <- 2
#' K <- 2
#' B <- 5
#' mean_dist <- 4
#' batch_dist <- 0.3
#' group_means <- seq(1, K) * mean_dist
#' batch_shift <- rnorm(B, mean = batch_dist, sd = batch_dist)
#' group_std_devs <- rep(2, K)
#' batch_scale <- rep(1.2, B)
#' group_weights <- rep(1 / K, K)
#' batch_weights <- rep(1 / B, B)
#' type <- "MVT"
#' group_dfs <- c(4, 7)
#' frac_known <- 0.3
#' permute_variables <- TRUE
#' scale_data <- FALSE
#' 
#' checkDataGenerationInputs(N,
#'   P,
#'   group_means,
#'   group_std_devs,
#'   batch_shift,
#'   batch_scale,
#'   group_weights,
#'   batch_weights,
#'   type,
#'   group_dfs,
#'   frac_known,
#'   permute_variables,
#'   scale_data
#' )
checkDataGenerationInputs <- function(N,
                                      P,
                                      group_means,
                                      group_std_devs,
                                      batch_shift,
                                      batch_scale,
                                      group_weights,
                                      batch_weights,
                                      type,
                                      group_dfs,
                                      frac_known,
                                      permute_variables,
                                      scale_data) {
  
  N_not_positive <- (N <= 0)
  N_not_whole_number <- (N %% 1) != 0
  if(N_not_positive || N_not_whole_number) {
    stop("N must be a positive whole number.")
  }
  
  P_not_positive <- (P <= 0)
  P_not_whole_number <- (P %% 1) != 0
  if(P_not_positive || P_not_whole_number) {
    stop("P must be a positive whole number.")
  }
  
  valid_density_passed <- (type == "MVN" || type == "MVT")
  if (!valid_density_passed) {
    stop("``type`` must be one of 'MVN' or 'MVT'.")
  }
  
  mvn_generated <- type == "MVN"
  mvt_generated <- type == "MVT"
  
  if (mvt_generated) {
    dfs_not_passed <- is.null(group_dfs)
    if (dfs_not_passed) {
      .err <- paste0(
        "If generating from a MVT, please pass values for the group degrees ",
        "of freedom (``group_dfs``)."
      )
      stop(.err)
    }
    dfs_too_small <- any(group_dfs < 1.0)
    if (dfs_too_small) {
      stop("Group degrees of freedom must be greater than or equal to 1.")
    }
  }
  
  # The number of batches to generate
  B <- length(batch_weights)
  
  wrong_length_batch_shift <- length(batch_shift) != B
  wrong_length_batch_scale <- length(batch_scale) != B
  
  # Allow for varying groups within batches and then record the number of groups
  # to generate
  varying_group_within_batch <- is.matrix(group_weights)
  if (varying_group_within_batch) {
    K <- nrow(group_weights)
    if (ncol(group_weights) != B) {
      stop("Number of columns in group weight matrix does not match the number of batches.")
    }
  } else {
    K <- length(group_weights)
  }
  
  wrong_length_group_means <- length(group_means) != K
  wrong_length_group_std_devs <- length(group_std_devs) != K
  
  
  if (any(
    wrong_length_group_means,
    wrong_length_group_std_devs,
    wrong_length_batch_shift,
    wrong_length_batch_scale
  )) {
    .err <- paste0(
      "Number of group and batch weights must match the number of ",
      "parameters in each of ``group_means``, ``group_std_devs``, ",
      "``batch_shift`` and ``batch_scale``."
    )
    stop(.err)
  }
  
  frac_known_not_in_0_1 <- (frac_known < 0.0 || frac_known > 1.0)
  if (frac_known_not_in_0_1) {
    stop("``frac_known`` must be a value between 0 and 1 (inclusive).")
  }
  
  batch_scale_too_small <- any(batch_scale < 1.0)
  if (batch_scale_too_small) {
    stop("All entries in ``batch_scale`` must be greater than or equal to 1.0.")
  }
  non_positive_std_devs <- any(group_std_devs <= 0.0)
  if (non_positive_std_devs) {
    stop("All entries of ``group_std_devs`` must be strictly positive.")
  }
  
  permute_variables_not_logical <- ! is.logical(permute_variables)
  if(permute_variables_not_logical) {
    stop("``permute_variable`` must be logical.")
  }
  
  scale_data_not_logical <- ! is.logical(scale_data)
  if(scale_data_not_logical) {
    stop("``scale_data`` must be logical.")
  }
}