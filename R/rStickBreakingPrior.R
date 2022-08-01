#!/usr/bin/Rscript
#' @title Random Draw From Stick Breaking Prior
#' @description Draw weights from the stick-breaking prior.
#' @param alpha The concentration parameter.
#' @param K The number of weights to generate.
#' @return A vector of component weights.
#' @export
#' @examples
#' weights <- rStickBreakingPrior(1, 50)
#' @importFrom stats rbeta
 rStickBreakingPrior <- function(alpha, K) {
  v <- stats::rbeta(K, 1, alpha)
  stick <- 1
  w <- rep(0, K)
  
  for (i in seq(1, K)) {
    w[i] <- v[i] * stick
    stick <- stick - w[i]
  }
  w
}
