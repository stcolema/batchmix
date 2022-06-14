#!/usr/bin/Rscript
#' @title Prepare initial values
#' @description Prepares user given values for input into the C++ function.
#' @param initial_class_means A $P x K$ matrix of initial values for the class
#' means. Defaults to draws from the prior distribution.
#' @param initial_class_covariance A $P x P x K$ array of initial values for the
#' class covariance matrices. Defaults to draws from the prior distribution.
#' @param initial_batch_shift A $P x B$ matrix of initial values for the batch
#' shift effect Defaults to draws from the prior distribution.
#' @param initial_batch_scale A $P x B$ matrix of initial values for the batch
#' scales Defaults to draws from the prior distribution.
#' @param initial_class_df A $K$ vector of initial values for the class degrees
#' of freedom. Defaults to draws from the prior distribution.
#' @param P Integer. The number of measurements for each sample in the dataset
#' being modelled.
#' @param K The number of classes/clusters being modelled.
#' @param B The number of batches being modelled.
#' @param type The type of mixture model used; one of "MVN" or "MVT".
#' @return A named list containing the different parameters.
prepareInitialParameters <- function(initial_class_means,
                                     initial_class_covariance,
                                     initial_batch_shift,
                                     initial_batch_scale,
                                     initial_class_df,
                                     P,
                                     K,
                                     B,
                                     type) {

  # Define placeholders. These will be used to pass to C++ if no initial values
  # are given. C++ isn't great for optional arguments, so this bypasses the
  # problem
  class_means <- matrix(-Inf, nrow = P, ncol = K)
  class_cov <- array(-Inf, dim = c(P, P, K))
  batch_shift <- matrix(-Inf, nrow = P, ncol = B)
  batch_scale <- matrix(-Inf, nrow = P, ncol = B)
  class_df <- rep(-Inf, K)

  class_mean_passed <- !is.null(initial_class_means)
  class_covariance_passed <- !is.null(initial_class_covariance)
  batch_shift_passed <- !is.null(initial_batch_shift)
  batch_scale_passed <- !is.null(initial_batch_scale)
  class_df_passed <- !is.null(initial_class_df)

  if (class_mean_passed) {
    not_a_matrix <- !is.matrix(initial_class_means)
    if (not_a_matrix) {
      stop("Initial class means must be either ``NULL`` or a P x K matrix.")
    }

    incorrect_dimensions <- ((nrow(initial_class_means) != P) &
      (ncol(initial_class_means) != K)
    )

    if (incorrect_dimensions) {
      stop("Initial class means must be either ``NULL`` or a P x K matrix.")
    }

    class_means <- initial_class_means
  }

  if (class_covariance_passed) {
    error_msg <- paste0(
      "Initial class covariances must be either ``NULL`` or ",
      "a P x P x K array"
    )

    not_an_array <- !is.array(initial_class_covariance)
    if (not_an_array) {
      stop(error_msg)
    }

    incorrect_dimensions <- any(dim(initial_class_covariance) != c(P, P, K))

    if (not_a_matrix | incorrect_dimensions) {
      stop(error_msg)
    }

    class_cov <- initial_class_covariance
  }

  if (class_df_passed) {
    wrong_type <- type != "MVT"
    if (wrong_type) {
      warning("Initial class degrees of freedom passed but model is MVN. Please check.")
    }

    error_msg <- paste0(
      "Initial class degrees of freedom must be either ",
      "``NULL`` or a K vector (such as ``c(3, 5)``."
    )

    not_a_vector <- !is.vector(initial_class_df)
    if (not_a_vector) {
      stop(error_msg)
    }

    incorrect_dimensions <- length(initial_class_df) != K

    if (incorrect_dimensions) {
      stop(error_msg)
    }

    class_df <- initial_class_df
  }


  if (batch_shift_passed) {
    error_msg <- paste0(
      "Initial batch shifts must be either ``NULL`` or a ",
      "P x B matrix."
    )

    not_a_matrix <- !is.matrix(initial_batch_shift)
    if (not_a_matrix) {
      stop(error_msg)
    }

    incorrect_dimensions <- ((nrow(initial_batch_shift) != P) &
      (ncol(initial_batch_shift) != B)
    )

    if (incorrect_dimensions) {
      stop(error_msg)
    }

    batch_shift <- initial_batch_shift
  }

  if (batch_scale_passed) {
    error_msg <- paste0(
      "Initial batch scales must be either ``NULL`` or a ",
      "P x B matrix."
    )

    not_a_matrix <- !is.matrix(initial_batch_scale)
    if (not_a_matrix) {
      stop(error_msg)
    }

    incorrect_dimensions <- ((nrow(initial_batch_scale) != P) &
      (ncol(initial_batch_scale) != B)
    )

    if (incorrect_dimensions) {
      stop(error_msg)
    }

    batch_scale <- initial_batch_scale
  }

  list(
    class_means = class_means,
    class_cov = class_cov,
    batch_shift = batch_shift,
    batch_scale = batch_scale,
    class_df = class_df
  )
}
