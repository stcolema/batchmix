#!/usr/bin/Rscript
#' @title Generate initial labels
#' @description For simulated data, generates an initial labelling for sampling.
#' @param alpha The mass in the stick breaking prior
#' @param K The number of classes available.
#' @param fixed The vector of 0s and 1s indicating which labels are to be held
#' fixed.
#' @param labels The initial labelling. Defaults to NULL.
#' @param X Optional N x P data matrix. If given and the model is fully
#' unsupervised (``fixed`` is all 0), the unfixed labels are seeded from a
#' k-means partition of ``X`` instead of a random draw from the stick-
#' breaking prior. A random draw is completely uninformed by the data, which
#' for well-separated clusters makes it easy for a chain to start (and,
#' since single-site Metropolis-Hastings allocation/mean updates essentially
#' never cross a low-density valley in a realistic run length, get
#' permanently stuck for the rest of the chain) in a poor local optimum -
#' this shows up as chains that individually look well-mixed but disagree
#' persistently with one another, i.e. a high rank-normalized Rhat that gets
#' \emph{worse}, not better, with a longer run (see
#' \code{\link{assessConvergence}}). Seeding from k-means gives every chain
#' a data-informed starting point instead, which is standard practice for
#' mixture-model MCMC initialisation (see e.g. Fruhwirth-Schnatter (2006),
#' "Finite Mixture and Markov Switching Models", Section 3.5) and does not
#' bias the posterior - it only changes where each chain starts, not its
#' stationary distribution. Falls back to the random draw if k-means fails
#' (e.g. fewer unfixed points than ``K``) or is not applicable (any fixed
#' points, or ``X`` not given).
#' @return An N vector of labels.
#' @export
#' @examples
#' N <- 500
#' P <- 2
#' K <- 2
#' B <- 5
#' mean_dist <- 4
#' batch_dist <- 0.3
#' cluster_means <- seq(1, K) * mean_dist
#' batch_shift <- rnorm(B, mean = batch_dist, sd = batch_dist)
#' std_dev <- rep(2, K)
#' batch_var <- rep(1.2, B)
#' cluster_weights <- rep(1 / K, K)
#' batch_weights <- rep(1 / B, B)
#'
#' my_data <- generateBatchData(
#'   N,
#'   P,
#'   cluster_means,
#'   std_dev,
#'   batch_shift,
#'   batch_var,
#'   cluster_weights,
#'   batch_weights
#' )
#'
#' initial_labels <- generateInitialLabels(1, K, my_data$fixed)
generateInitialLabels <- function(alpha, K, fixed, labels = NULL, X = NULL) {
  N <- length(fixed)
  if (is.null(labels)) {
    labels <- rep(0, N)
  }

  N_fixed <- sum(fixed)
  N_unfixed <- N - N_fixed
  labels_available <- seq(1, K)

  fully_unsupervised <- N_fixed == 0
  try_kmeans <- fully_unsupervised && !is.null(X) && N_unfixed >= K

  if (try_kmeans) {
    kmeans_labels <- tryCatch(
      stats::kmeans(X[fixed == 0, , drop = FALSE], centers = K, nstart = 10)$cluster,
      error = function(e) NULL
    )
    if (!is.null(kmeans_labels)) {
      labels[fixed == 0] <- kmeans_labels
      return(labels)
    }
  }

  # Draw prior weights from the stick breaking prior
  w <- rStickBreakingPrior(alpha, K)

  labels[fixed == 0] <- sample(labels_available, N_unfixed,
    replace = TRUE,
    prob = w
  )

  labels
}
