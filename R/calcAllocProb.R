#!/usr/bin/Rscript
#' @title Calculate allocation probabilities
#' @description Calculate the empirical allocation probability for each class
#' based on the sampled allocation probabilities.
#' @param mcmc_samples Output from ``batchSemiSupervisedMixtureModel``.
#' @param burn The number of samples to discard.
#' @param method The point estimate to use. ``method = 'mean'`` or
#' ``method = 'median'``. ``'median'`` is the default.
#' @return An N x K matrix of class probabilities.
#' @export
#' @examples
#' 
#' # Data in matrix format
#' X <- matrix(c(rnorm(100, 0, 1), rnorm(100, 3, 1)), ncol = 2, byrow = TRUE)
#' 
#' # Initial labelling
#' labels <- c(rep(1, 10), 
#'   sample(c(1,2), size = 40, replace = TRUE), 
#'   rep(2, 10), 
#'   sample(c(1,2), size = 40, replace = TRUE)
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
#'
#' # MCMC samples and BIC vector
#' samples <- batchSemiSupervisedMixtureModel(X, R, thin, labels, fixed, batch_vec, "MVN")
#' 
#' # Burn in
#' burn <- 20
#' eff_burn <- burn / thin
#' 
#' # Probability across classes
#' probs <- calcAllocProb(samples, burn = burn)
#' 
calcAllocProb <- function(mcmc_samples, burn = 0, method = "median") {
  R <- mcmc_samples$R
  thin <- mcmc_samples$thin
  .alloc <- mcmc_samples$alloc
  
  if(burn > 0) {
    if(burn > R) {
      stop("Burn in exceeds number of iterations run.")
    }
    
    eff_burn <- floor(burn / thin)
    dropped_samples <- seq(1, eff_burn)
    .alloc <- .alloc[, , -dropped_samples]
  }
  
  probs <- NULL
  
  if(method == "median") {
    probs <- apply(.alloc, c(1, 2), median)
  }
  if(method == "mean") {
    probs <- rowSums(.alloc, dims = 2) / dim(.alloc)[3]
  }
  if(length(probs) == 1) {
    if(is.null(probs)) {
      stop("``method`` must be one of 'mean' or 'median'")
    }
  }
  probs
}