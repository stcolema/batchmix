#!/usr/bin/Rscript
#' @title Plot likelihoods
#' @description Plots the model fit for multiple chains.
#' @param mcmc_outputs The output from ``runMCMCChains``.
#' @param choice The model fit score to use. Must be one of
#' ``'observed_likelihood'``, ``'complete_likelihood'`` or ``'BIC'``. Defaults
#' to ``'complete_likelihood'``.
#' @param colour_by_chain Logical indcating if plots should be coloured by chain
#' or all the same colour. Defaults to ``TRUE``.
#' @return A ggplot2 object. Line plot of likelihood across iteration.
#' @importFrom ggplot2 ggplot aes_string geom_line
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
#' R <- 1000
#' thin <- 50
#' n_chains <- 4
#'
#' # MCMC samples and BIC vector
#' samples <- runMCMCChains(X, n_chains, R, thin, labels, fixed, batch_vec, "MVN")
#'
#' p <- plotLikelihoods(samples)
#' 
plotLikelihoods <- function(mcmc_outputs, 
                            choice = "complete_likelihood",
                            colour_by_chain = TRUE
                            ) {
  lkl_lst <- lapply(mcmc_outputs, getLikelihood, choice)

  n_chains <- length(lkl_lst)
  for(ii in seq(1, n_chains)) {
    lkl_lst[[ii]]$Chain <- mcmc_outputs[[ii]]$Chain
  }

  lkl_df <- do.call(rbind, lkl_lst)
  lkl_df$Chain <- factor(lkl_df$Chain)

  if(colour_by_chain) {
    
    p <- ggplot2::ggplot(
      data = lkl_df,
      mapping = ggplot2::aes_string(
        x = "iteration",
        y = choice,
        colour = "Chain"
      )
    ) +
      ggplot2::geom_line()
  } else {
  p <- ggplot2::ggplot(
    data = lkl_df,
    mapping = ggplot2::aes_string(
      x = "iteration",
      y = choice,
      group = "Chain"
    )
  ) +
    ggplot2::geom_line()
}
  p
}
