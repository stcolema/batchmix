#' @title Bayesian Mixture Modelling for Joint Model-Based
#' Clustering/Classification and Batch Correction
#' @description Semi-supervised and unsupervised Bayesian mixture models that
#' simultaneously infer the cluster/class structure and a batch correction.
#' Densities available are the multivariate normal and the multivariate t.
#' The model sampler is implemented in C++. This package is aimed at analysis of
#' low-dimensional data generated across several batches. See
#' \href{https://doi.org/10.1101/2022.01.14.476352}{Coleman et al. (2022)} for
#' details of the model.
#'
#' @section Defaults that changed in this release:
#' Partial pooling is on by default: each batch has its own mixture weights
#' (\code{batch_weight_prior = "partial_pooling"}, or \code{"global"} for one
#' batch) and the concentration of the batch-scale prior is estimated
#' (\code{sample_s_scale = TRUE}). For the earlier fixed-hyperparameter model
#' use \code{batch_weight_prior = "global", sample_s_scale = FALSE}. Under
#' partial pooling \code{alpha}/\code{concentration} are unused, so choose
#' \code{K_max} near the number of clusters expected. The batch-corrected
#' data are now the posterior mean of each item's batch-free signal, several
#' sampler corrections change fixed-seed output, models should be compared
#' with \code{\link{calcBICM}}, and \code{salso} is no longer needed (see
#' \code{\link{minVI}}). \code{news(package = "batchmix")} has the details.
#'
#' @section Where to start:
#' \code{\link{fitBatchMix}} is the main entry point. Read the vignettes
#' (\code{browseVignettes("batchmix")}): \code{batchmix_workflow} (start
#' here), \code{batch_weight_priors}, \code{covariance_models} and
#' \code{probit_missing_censored}.
#' @name batchmix-package
#' @aliases batchmix
#' @docType package
#' @author Stephen Coleman <stcolema@tcd.ie>, Paul D.W. Kirk, Chris Wallace
#' @keywords package
#'
#' @importFrom ggplot2 aes facet_grid facet_wrap geom_boxplot geom_line geom_point
#' ggplot label_both labeller labs
#' @importFrom stats cutree as.dist hclust median rbeta rchisq rnorm
#' @importFrom tidyr contains pivot_longer
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib batchmix
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
#' n_iter <- 1000
#' thin <- 50
#'
#' # Classification
#' samples <- runBatchMix(X,
#'   n_iter,
#'   thin,
#'   batch_vec,
#'   "MVN",
#'   initial_labels = labels,
#'   fixed = fixed,
#' )
#'
#' # Clustering
#' samples <- runBatchMix(X, n_iter, thin, batch_vec, "MVT")
#'
NULL
