#!/usr/bin/Rscript
#' @title Plot sampled vector parameter
#' @description Plot the sampled values for a sampled vector from the output of
#' the ``mixtureModel`` function. Not recommended for large B or P.
#' @param samples The output of the ``mixtureModel`` function.
#' @param parameter The name of the parameter to be plotted (a string).
#' @param R The number of iterations run. Defaults to the number of samples for
#' the cluster membership.
#' @param thin The thinning factor of the sampler. Defaults to 1. 
#' @param burn_in The samples at the beginning of the chain to drop. Defaults to 0.
#' @return A ggplot object of the values in each sampled batch mean per iteration.
#' @export
#' @examples
#' # Data in matrix format
#' X <- matrix(c(rnorm(100, 0, 1), rnorm(100, 3, 1)), ncol = 2, byrow = TRUE)
#' 
#' # Observed batches represented by integers
#' batch_vec <- sample(1:5, size = 100, replace = TRUE)
#' 
#' # MCMC iterations (this is too low for real use)
#' R <- 50
#' thin <- 1
#'
#' # MCMC samples and BIC vector
#' samples <- batchUnsupervisedMixtureModel(X, R, thin, batch_vec, "MVN")
#' 
#' # Plot the sampled value of the cluster means against MCMC iteration 
#' parameter <- "means"
#' plotSampledParameter(samples, parameter, R, thin)
#' @importFrom ggplot2 ggplot aes_string geom_point facet_wrap labs
plotSampledParameter <- function(samples, parameter, R = NULL, thin = 1, burn_in = 0) {
  n_param <- dim(samples[[parameter]])[2]
  P <- dim(samples[[parameter]])[1]
  
  if (is.null(R)) {
    R <- nrow(samples$samples)
  }
  
  # Check that the values of R and thin make sense
  if(floor(R/thin) != nrow(samples$samples)){
    stop("The ratio of R to thin does not match the number of samples present.")
  }
  
  sampled_parameter <- getSampledBatchShift(samples[[parameter]], n_param, P, R = R, thin = thin)
  
  sampled_parameter <- sampled_parameter[sampled_parameter$Iteration > burn_in, ]
  
  p <- ggplot2::ggplot(sampled_parameter, 
      ggplot2::aes_string(x = "Iteration", y = "value")
    ) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~name, ncol = P) +
    ggplot2::labs(
      title = parameter,
      x = "MCMC iteration",
      y = "Sampled value"
    )
  
  p
}
