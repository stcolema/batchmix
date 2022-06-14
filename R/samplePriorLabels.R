#!/usr/bin/Rscript
#' @title Sample prior labels
#' @description Generate labels from the stick-breaking prior.
#' @param alpha The concentration parameter for the stick-breaking prior.
#' @param K The number of components to include (the upper bound on the number of unique labels generated).
#' @param N The number of labels to generate.
#' @return A vector of labels.
#' @export
#' @examples
#' initial_labels <- samplePriorLabels(1, 50, 100)
samplePriorLabels <- function(alpha, K, N) {
  w <- stickBreakingPrior(alpha, K)
  initial_labels <- sample(seq(1, K), N, replace = TRUE, prob = w)
}
