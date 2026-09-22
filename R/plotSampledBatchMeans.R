#!/usr/bin/Rscript
#' @title Plot sampled batch means
#' @description Plot the sampled values for the batch mean shifts in each
#' dimension from the output of the mixture model functions. Not recommended
#' for large B or P.
#' @param samples The output of the ``batchUnsupervisedMixtureModel`` or
#' ``batchSemiSupervisedMixtureModel`` functions.
#' @param burn_in The samples at the beginning of the chain to drop. Defaults to 0.
#' @return A ggplot object of the values in each sampled batch mean per iteration.
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
#' # Plot the sampled value of the batch mean shift against MCMC iteration
#' plotSampledBatchMeans(samples)
#'
#' @importFrom ggplot2 ggplot aes geom_point facet_grid labs labeller label_both
plotSampledBatchMeans <- function(samples, burn_in = 0) {
  B <- samples$B
  P <- samples$P

  n_iter <- samples$n_iter
  thin <- samples$thin

  # Check that the values of n_iter and thin make sense
  if (floor(n_iter / thin) != nrow(samples$samples)) {
    stop("The ratio of n_iter to thin does not match the number of samples present.")
  }

  sampled_batch_shift <- getSampledBatchShift(samples$batch_shift, B, P,
    n_iter = n_iter,
    thin = thin
  )

  # Remove the warm-up samples
  sampled_batch_shift <- sampled_batch_shift[
    sampled_batch_shift$Iteration > burn_in,
  ]

  # Make a ggplot2 object
  p <- ggplot2::ggplot(
    sampled_batch_shift,
    ggplot2::aes(x = Iteration, y = value)
  ) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(Batch ~ Dimension,
      labeller = ggplot2::labeller(
        Batch = ggplot2::label_both,
        Dimension = ggplot2::label_both
      )
    ) +
    ggplot2::labs(
      title = "Batch mean shift",
      x = "MCMC iteration",
      y = "Sampled value"
    )

  p
}
