#' @title Plot a prior or posterior predictive check
#' @description Visualises replicated datasets (from
#' \code{\link{simulatePriorPredictive}} or
#' \code{\link{simulatePosteriorPredictive}}) against the observed data, in
#' either of the two standard forms used in a Bayesian workflow (Gelman,
#' Carlin, Stern, Dunson, Vehtari & Rubin, 2013, \emph{Bayesian Data
#' Analysis}, 3rd ed., Section 6.3-6.4; Gabry, Simpson, Vehtari, Betancourt &
#' Gelman, 2019, "Visualization in Bayesian Workflow", JRSS-A 182(2)):
#' \itemize{
#'   \item \code{style = "density"}: overlays the marginal density of one
#'   observed column against the same column's density in each replicate
#'   dataset (analogous to \code{bayesplot::ppc_dens_overlay()}) - a
#'   graphical check of overall distributional shape.
#'   \item \code{style = "statistic"}: computes a scalar test statistic (e.g.
#'   the mean or SD of a column, optionally within one batch) on the
#'   observed data and on every replicate, and histograms the replicate
#'   values with the observed value marked - the classical (BDA3 Eq. 6.9)
#'   posterior/prior predictive check, whose tail-area gives a Bayesian
#'   p-value. A value near 0 or 1 flags a real discrepancy the model fails
#'   to capture (e.g. a batch effect the model is not removing); values
#'   scattered through the middle of the range are the "nothing detected"
#'   outcome, as expected for most statistics most of the time.
#' }
#' @param X The observed N x P data matrix.
#' @param replicates A list as returned by \code{\link{simulatePriorPredictive}}
#' or \code{\link{simulatePosteriorPredictive}} (each element has an
#' \code{X} entry).
#' @param style One of \code{"density"} or \code{"statistic"}. For a binary/
#' probit column (\code{type = "MVN_MIXED"}, \code{column_type == 1} in the
#' original fit), use \code{"statistic"} with \code{statistic = mean} (the
#' proportion of 1s) rather than \code{"density"}, which assumes a
#' continuous column - see the "binary column" warning below.
#' @param column Which column of X to check (used by both styles).
#' @param censor_code Optional N x P matrix, non-zero where \code{X} was
#' censored (see \code{\link{batchSemiSupervisedMixtureModel}}). When
#' supplied, censored cells in \code{column} are excluded from the
#' *observed* side of the comparison: \code{X}'s recorded value there is a
#' censoring bound, not the true value, so it is not a fair comparison
#' target for a replicate's freely-simulated draw (every replicate is fully
#' observed - see \code{\link{simulatePriorPredictive}}/
#' \code{\link{simulatePosteriorPredictive}}). \code{NA} entries in
#' \code{X} (ordinary missingness) are always excluded from the observed
#' side this way too, \code{censor_code} or not.
#' @param statistic Only used if \code{style = "statistic"}: a function
#' taking a numeric vector and returning a single number (e.g. \code{mean},
#' \code{stats::sd}, or a closure over a batch-specific subset - see
#' examples). Applied to the non-\code{NA} entries only.
#' @param statistic_name Only used if \code{style = "statistic"}: a label
#' for the x-axis.
#' @param max_replicates Only used if \code{style = "density"}: caps how many
#' replicate densities are drawn, to keep the plot legible.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' \donttest{
#' X <- matrix(c(rnorm(100, 0, 1), rnorm(100, 3, 1)), ncol = 2, byrow = TRUE)
#' batch_vec <- sample(seq(1, 3), replace = TRUE, size = 100)
#'
#' prior_sims <- simulatePriorPredictive(X, batch_vec, K = 2, type = "MVN", n_datasets = 20)
#' plotPredictiveCheck(X, prior_sims, style = "density", column = 1)
#' plotPredictiveCheck(X, prior_sims,
#'   style = "statistic", column = 1,
#'   statistic = stats::sd, statistic_name = "SD(column 1)"
#' )
#' }
plotPredictiveCheck <- function(X,
                                replicates,
                                style = c("density", "statistic"),
                                column = 1,
                                censor_code = NULL,
                                statistic = mean,
                                statistic_name = "statistic",
                                max_replicates = 50) {
  style <- match.arg(style)

  # The observed side of the comparison: X's own NAs (ordinary missingness)
  # are always excluded; censor_code additionally excludes cells whose
  # recorded value is only a censoring bound, not the true value - see the
  # censor_code documentation above.
  obs_values <- X[, column]
  if (!is.null(censor_code)) {
    obs_values[censor_code[, column] != 0] <- NA_real_
  }

  is_binary_like <- all(obs_values %in% c(0, 1, NA))
  if (style == "density" && is_binary_like) {
    warning(
      "Column ", column, " looks binary (every observed value is 0, 1 or NA). ",
      "geom_density() assumes a continuous column; style = 'statistic' with ",
      "statistic = mean (the proportion of 1s) is the appropriate check here."
    )
  }

  if (style == "density") {
    replicates <- replicates[seq_len(min(length(replicates), max_replicates))]

    obs_df <- data.frame(value = obs_values)
    rep_df <- do.call(rbind, lapply(seq_along(replicates), function(i) {
      data.frame(value = replicates[[i]]$X[, column], draw = i)
    }))

    p <- ggplot2::ggplot() +
      ggplot2::geom_density(
        data = rep_df,
        ggplot2::aes(x = .data$value, group = .data$draw),
        colour = "steelblue", alpha = 0.25, linewidth = 0.3
      ) +
      ggplot2::geom_density(
        data = obs_df, ggplot2::aes(x = .data$value),
        colour = "black", linewidth = 1
      ) +
      ggplot2::labs(
        title = paste0("Predictive check: column ", column),
        subtitle = paste0(length(replicates), " replicates (blue) vs. observed (black)"),
        x = paste0("Column ", column), y = "Density"
      ) +
      ggplot2::theme_minimal()

    return(p)
  }

  obs_stat <- statistic(obs_values[!is.na(obs_values)])
  rep_stats <- vapply(replicates, function(r) {
    rep_col <- r$X[, column]
    statistic(rep_col[!is.na(rep_col)])
  }, numeric(1))

  p_value <- mean(rep_stats >= obs_stat)

  ggplot2::ggplot(data.frame(value = rep_stats), ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = obs_stat, colour = "black", linewidth = 1) +
    ggplot2::labs(
      title = paste0("Predictive check: ", statistic_name, ", column ", column),
      subtitle = paste0(
        "Observed = ", signif(obs_stat, 4),
        " (black line); Bayesian p-value = P(replicate >= observed) = ",
        signif(p_value, 3)
      ),
      x = statistic_name, y = "Count across replicates"
    ) +
    ggplot2::theme_minimal()
}
