#' @title Bayesian Mixture Modelling for Joint Model-Based
#' Clustering/Classification and Batch Correction
#' @description Semi-supervised and unsupervised Bayesian mixture models that
#' simultaneously infer the cluster/class structure and a batch correction.
#' Densities available are the multivariate normal and the multivariate t.
#' The model sampler is implemented in C++. This package is aimed at analysis of
#' low-dimensional data generated across several batches. See
#' (Coleman et al. (2022))[https://doi.org/10.1101/2022.01.14.476352] for
#' details of the model.
#' @name batchmix-package
#' @aliases batchmix, BatchMixtureModel
#' @docType package
#' @author Stephen Coleman <stcolema@tcd.ie>, Paul D.W. Kirk, Chris Wallace
#' @keywords package
#'
#' @importFrom ggplot2 aes_string facet_grid facet_wrap geom_line geom_point
#' ggplot label_both labeller labs
#' @importFrom stats median rbeta rchisq rnorm
#' @importFrom tidyr contains pivot_longer
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib BatchMixtureModel
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
#' # Which labels are observed
#' fixed <- c(rep(1, 10), rep(0, 40), rep(1, 10), rep(0, 40))
#'
#' # Batch
#' batch_vec <- sample(seq(1, 5), replace = TRUE, size = 100)
#'
#' # Sampling parameters
#' R <- 1000
#' thin <- 50
#'
#' # Classification
#' samples <- runBatchMix(X,
#'   R,
#'   thin,
#'   batch_vec,
#'   "MVN",
#'   initial_labels = labels,
#'   fixed = fixed,
#' )
#'
#' # Clustering
#' samples <- runBatchMix(X, R, thin, batch_vec, "MVT")
#'
NULL
