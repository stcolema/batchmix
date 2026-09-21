#' @title Minimum VI lower bound
#' @description Local implementation of S. Wade's `minVI` function from their
#' `mcclust.ext` package (available from github).
#' Reimplemented here to avoid dependency on a non-CRAN package. Computes the
#' lower bound to the posterior expected Variation of Information. For full
#' details please see the aforementioned package and Wade and Ghahramani, 2018,
#' 'Bayesian Cluster Analysis: Point Estimation and Credible Balls (with
#' Discussion)'.
#' @param cls A clustering for which the lower bound of the Variation of
#' Information is calculated.
#' @param psm The posterior similarity matrix which `cls` is a summary thereof.
#' @returns A vector of the the lower bound of the Variation of Information
#' for
#' @examples
#' \dontrun{
#' # MCMC samples and BIC vector
#' mcmc_outputs <- runMCMCChains(
#'   X,
#'   n_chains,
#'   R,
#'   thin,
#'   batch_vec,
#'   type
#' )
#'
#' # Note that in this toy example we have not applied a burn in
#' psm <- createSimilarityMat(mcmc_outputs[[1]]$samples)
#' VI.lb(mcmc_outputs[[1]]$samples[1, ], psm)
#' }
VI.lb <- function(cls, psm) {
  if (is.vector(cls)) {
    cls <- t(cls)
  }
  n <- nrow(psm)
  n_inds <- seq(1, n)

  VI.lb.compute <- function(c) {
    f <- 0
    n_inds <- seq(1, n)
    for (i in n_inds) {
      ind <- (c == c[i])
      f <- f + (log2(sum(ind))
      + log2(sum(psm[i, ]))
        - 2 * log2(sum(ind * psm[i, ]))
      ) / n
    }
    return(f)
  }
  output <- apply(cls, 1, VI.lb.compute)
  output
}
