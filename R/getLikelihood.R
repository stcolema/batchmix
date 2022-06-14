#!/usr/bin/Rscript
#' @title Get likelihood
#' @description Extracts the model fit score from the mixture model output.
#' @param mcmc_output The output from the mixture model.
#' @param choice The model fit score to use. Must be one of
#' ``'observed_likelihood'``, ``'complete_likelihood'`` or ``'BIC'``. Defaults
#' to ``'complete_likelihood'``.
#' @return A data.frame containing the model fit score of choice and the
#' iteration.
#' @export
#' @examples
#' 
#' # Data in a matrix format
#' X <- matrix(c(rnorm(100, 0, 1), rnorm(100, 3, 1)), ncol = 2, byrow = TRUE)
#'
#' # Batch
#' batch_vec <- sample(seq(1, 5), replace = TRUE, size = 100)
#'
#' # Sampling parameters
#' R <- 100
#' thin <- 5
#'
#' # MCMC samples and BIC vector
#' samples <- batchUnsupervisedMixtureModel(X, R, thin, batch_vec, "MVN")
#'
#' lkl_df <- getLikelihood(samples)
#' 
getLikelihood <- function(mcmc_output, choice = "complete_likelihood") {
  R <- mcmc_output$R
  thin <- mcmc_output$thin
  burn <- mcmc_output$burn
  first_recorded_iter <- burn + thin

  iters <- seq(first_recorded_iter, R, by = thin)

  invalid_choice <- !(
    choice %in% c(
      "observed_likelihood",
      "complete_likelihood",
      "BIC"
    )
  )

  if (invalid_choice) {
    stop("Choice must be one of 'observed_likelihood', 'complete_likelihood' or 'BIC'.")
  }

  sampled_likelihoods <- mcmc_output[[choice]]

  lkl_df <- data.frame(score = sampled_likelihoods, iteration = iters)
  colnames(lkl_df) <- c(choice, "iteration")

  lkl_df
}
