#!/usr/bin/Rscript
#' @title Random Draw From Stick Breaking Prior
#' @description Draw weights from the stick-breaking prior of a Dirichlet process
#' with concentration \code{alpha}, truncated at \code{K} components (the last
#' weight takes the remaining mass), so that the weights sum to one.
#' @param alpha The concentration parameter.
#' @param K The number of weights to generate.
#' @return A vector of \code{K} component weights summing to one.
#' @export
#' @examples
#' weights <- rStickBreakingPrior(1, 50)
#' @importFrom stats rbeta
rStickBreakingPrior <- function(alpha, K) {
  v <- stats::rbeta(K, 1, alpha)
  # Truncation at K components (Ishwaran & James, 2001, JASA 96(453)): the last
  # break takes whatever remains of the stick, so that the weights sum to one
  v[K] <- 1
  stick <- 1
  w <- rep(0, K)

  for (i in seq(1, K)) {
    w[i] <- v[i] * stick
    stick <- stick - w[i]
  }
  w
}
