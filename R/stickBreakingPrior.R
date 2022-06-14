#!/usr/bin/Rscript
#' @title Stick breaking prior
#' @description Draw weights from the stick-breaking prior.
#' @param alpha The concentration parameter.
#' @param K The number of weights to generate.
#' @return A vector of weights.
#' @export
#' @examples
#' weights <- stickBreakingPrior(1, 50)
#' @importFrom stats rbeta
stickBreakingPrior <- function(alpha, K) {
  v <- stats::rbeta(K, alpha, 1)
  stick <- 1
  w <- rep(0, K)
  
  for (i in seq(1, K)) {
    w[i] <- v[i] * stick
    stick <- stick - w[i]
  }
  w
}
