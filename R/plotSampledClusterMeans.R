#!/usr/bin/Rscript
#' @title Plot sampled cluster means
#' @description Plot the sampled values for the cluster means in each
#' dimension from the output of the mixture model functions. Not recommended
#' for large K or P.
#' @param samples The output of the ``batchUnsupervisedMixtureModel`` or
#' ``batchSemiSupervisedMixtureModel`` functions.
#' @param burn_in The samples at the beginning of the chain to drop. Defaults to 0.
#' @return A ggplot object of the values in each sampled cluster mean per iteration.
#' @export
#' @examples
#'
#' # Data in matrix format
#' X <- matrix(c(rnorm(100, 0, 1), rnorm(100, 3, 1)), ncol = 2, byrow = TRUE)
#'
#' # Observed batches represented by integers
#' batch_vec <- sample(seq(1, 5), size = 100, replace = TRUE)
#'
#' # MCMC iterations (this is too low for real use)
#' n_iter <- 100
#' thin <- 5
#'
#' # MCMC samples and BIC vector
#' samples <- runBatchMix(X, n_iter, thin, batch_vec, "MVN")
#'
#' # Plot the sampled value of the cluster means against MCMC iteration
#' plotSampledClusterMeans(samples)
#'
#' @importFrom ggplot2 ggplot aes geom_point facet_grid labs labeller label_both
plotSampledClusterMeans <- function(samples, burn_in = 0) {
  K <- samples$K_max
  P <- samples$P

  n_iter <- samples$n_iter
  thin <- samples$thin

  # Check that the values of n_iter and thin make sense
  if (floor(n_iter / thin) != nrow(samples$samples)) {
    stop("The ratio of n_iter to thin does not match the number of samples present.")
  }

  sampled_cluster_means <- getSampledClusterMeans(samples$means, K, P,
    n_iter = n_iter,
    thin = thin
  )

  # Remove warm-up samples
  sampled_cluster_means <- sampled_cluster_means[
    sampled_cluster_means$Iteration > burn_in,
  ]

  # Create the plot as aggplot2 object
  p <- ggplot2::ggplot(
    sampled_cluster_means,
    ggplot2::aes(x = Iteration, y = value)
  ) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(Cluster ~ Dimension,
      labeller = ggplot2::labeller(
        Cluster = ggplot2::label_both,
        Dimension = ggplot2::label_both
      )
    ) +
    ggplot2::labs(
      title = "Cluster means",
      x = "MCMC iteration",
      y = "Sampled value"
    )

  p
}
