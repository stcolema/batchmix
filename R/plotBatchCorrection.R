#' @title Visualise the effect of batch correction
#' @description The natural visual complement to the numerical batch-effect
#' checks (\code{\link{plotPredictiveCheck}}, the per-batch statistic
#' example in the \code{bayesian_workflow} vignette): a before/after scatter
#' (for \code{P >= 2}, using \code{columns}) or density (for \code{P == 1})
#' plot of the raw data against the point-estimate batch-corrected data,
#' coloured by batch. A correction that is working should visibly pull the
#' batches together in the "after" panel while leaving the cluster
#' structure (shape/relative position) recognisable - if it does not, or if
#' it collapses genuine cluster differences along with the batch
#' differences, that is worth investigating further (e.g. with
#' \code{\link{plotPredictiveCheck}}'s per-batch statistic check).
#' @param mcmc_output Output of \code{\link{batchSemiSupervisedMixtureModel}}/
#' \code{\link{runBatchMix}} (\code{type} one of 'MVN', 'MVT', 'MVN_LKJ' -
#' 'MVN_MIXED' is not currently supported here).
#' @param X The data matrix used to fit \code{mcmc_output}.
#' @param batch_vec The batch label used to fit \code{mcmc_output}.
#' @param burn Number of iterations to discard as burn-in before computing
#' the point-estimate correction (see \code{\link{processMCMCChain}}).
#' @param point_estimate_method Passed to \code{\link{processMCMCChain}}.
#' @param columns Only used if \code{ncol(X) >= 2}: which two columns to
#' plot.
#' @param show_raw Logical; if \code{TRUE} (the default), the plot is
#' faceted into "Raw" and "Batch-corrected" panels, as above. If
#' \code{FALSE}, only the batch-corrected panel is drawn (no facet, no
#' "Raw" panel) - useful once the raw-vs-corrected comparison has already
#' been made and only the corrected data is of interest, e.g. for a report
#' figure or downstream visual inspection of the data actually used for
#' clustering.
#' @return A \code{ggplot} object (for \code{P >= 2}, a scatter plot; for
#' \code{P == 1}, a density plot; faceted into "Raw"/"Batch-corrected"
#' panels unless \code{show_raw = FALSE}, in which case only the
#' batch-corrected panel is returned).
#' @export
#' @examples
#' \donttest{
#' X <- matrix(c(rnorm(100, 0, 1), rnorm(100, 3, 1)), ncol = 2, byrow = TRUE)
#' batch_vec <- sample(seq(1, 3), replace = TRUE, size = 100)
#' mcmc_out <- runBatchMix(X, 1000, 50, batch_vec, "MVN", verbose = FALSE)
#' plotBatchCorrection(mcmc_out, X, batch_vec, burn = 250)
#' plotBatchCorrection(mcmc_out, X, batch_vec, burn = 250, show_raw = FALSE)
#' }
plotBatchCorrection <- function(mcmc_output,
                                X,
                                batch_vec,
                                burn = 0,
                                point_estimate_method = "median",
                                columns = c(1, 2),
                                show_raw = TRUE) {
  if (isTRUE(mcmc_output$type == "MVN_MIXED")) {
    stop("plotBatchCorrection() does not support type = 'MVN_MIXED' yet.")
  }

  processed <- processMCMCChain(mcmc_output, burn, point_estimate_method)
  X_corrected <- processed$inferred_dataset

  P <- ncol(X)
  batch_fac <- factor(batch_vec)

  if (P == 1) {
    df <- data.frame(value = X_corrected[, 1], batch = batch_fac, panel = "Batch-corrected")
    if (show_raw) {
      df <- rbind(data.frame(value = X[, 1], batch = batch_fac, panel = "Raw"), df)
      df$panel <- factor(df$panel, levels = c("Raw", "Batch-corrected"))
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, colour = .data$batch)) +
      ggplot2::geom_density(linewidth = 0.8) +
      ggplot2::labs(title = "Effect of batch correction", x = "Column 1", colour = "Batch") +
      ggplot2::theme_minimal()

    if (show_raw) p <- p + ggplot2::facet_wrap(~panel)

    return(p)
  }

  col1 <- columns[1]
  col2 <- columns[2]

  df <- data.frame(x = X_corrected[, col1], y = X_corrected[, col2], batch = batch_fac, panel = "Batch-corrected")
  if (show_raw) {
    df <- rbind(data.frame(x = X[, col1], y = X[, col2], batch = batch_fac, panel = "Raw"), df)
    df$panel <- factor(df$panel, levels = c("Raw", "Batch-corrected"))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$batch)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::labs(
      title = "Effect of batch correction",
      x = paste0("Column ", col1), y = paste0("Column ", col2), colour = "Batch"
    ) +
    ggplot2::theme_minimal()

  if (show_raw) p <- p + ggplot2::facet_wrap(~panel)

  p
}
