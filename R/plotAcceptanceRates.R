#!/usr/bin/Rscript
#' @title Plot acceptance rates
#' @description Plot the acceptance rates for the parameters sampled in a 
#' Metropolis-Hastings step. Aiming for acceptance rates in [0.2, 0.5] for the 
#' class means, the batch effect on location and scale is a good rule of thumb. 
#' The class covariance should be in [0.35, 0.8] based on the authors' 
#' experience. The class degree of freedom appears to be prone to high 
#' acceptance rates, but aim to keep this above 0.2 at a minimum.
#' @param mcmc_lst The output of the ``runMCMCChains`` function.
#' @return A ggplot object of the boxplots of acceptance rates for each 
#' parameter across chains.
#' @export
#' @examples
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
#' R <- 500
#' thin <- 10
#' n_chains <- 4
#'
#' # MCMC samples and BIC vector
#' mcmc_lst <- runMCMCChains(X, n_chains, R, thin, batch_vec, "MVN", 
#'   initial_labels = labels, 
#'   fixed = fixed
#' )
#' 
#' # Plot the acceptance rate of each parameter in the 4 chains
#' plotAcceptanceRates(mcmc_lst)
#' 
#' @importFrom ggplot2 ggplot aes_string geom_boxplot
#' @importFrom tidyr pivot_longer any_of
plotAcceptanceRates <- function(mcmc_lst) {
  type <- mcmc_lst[[1]]$type
  
  acceptance_rates <- do.call(rbind,
    lapply(mcmc_lst, collectAcceptanceRates, type)
  )
  
  acceptance_rates$Chain <- seq(1, length(mcmc_lst))
  
  plot_df <- tidyr::pivot_longer(acceptance_rates, -tidyr::any_of("Chain"), 
    names_to = "Parameter",
    values_to = "Acceptance_rate"
  ) 
  
  p_out <- ggplot2::ggplot(plot_df, 
    ggplot2::aes_string(x = "Parameter", y = "Acceptance_rate")
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::ylim(0, 1)
  
  p_out
}