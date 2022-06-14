#!/usr/bin/Rscript
#' @title Predict class
#' @description Predicts a final class for each item given a matrix of 
#' allocation probabilities.
#' @param prob Output from the ``calcAllocProb`` function, a N x K matrix of 
#' allocation probabilities.
#' @return An N vector of class allocations.
#' @export
#' @examples
#' 
#' # Data in a matrix format
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
#' samples <- batchSemiSupervisedMixtureModel(X, 
#'   R, 
#'   thin, 
#'   labels, 
#'   fixed, 
#'   batch_vec, 
#'   "MVN"
#' )
#' 
#' # Burn in
#' burn <- 200
#' eff_burn <- burn / thin 
#' 
#' # Probability across classes
#' probs <- calcAllocProb(samples, burn = burn)
#' 
#' # Predict the class
#' preds <- predictClass(probs)
#' 
predictClass <- function(prob) {
  
  pred_cl <- apply(prob, 1, which.max)
  
  pred_cl
}