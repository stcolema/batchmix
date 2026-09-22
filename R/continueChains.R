#' @title Continue chains
#' @description Continues sampling from a list of previous chains.
#' @param mcmc_output Chains to be continued.
#' @param X Data to cluster as a matrix with the items to cluster held in rows.
#' @param fixed The indicator vector for which labels are observed.
#' @param batch_vec The vector of the batch labels for the data.
#' @param n_iter The number of iterations to run in this continuation (thinning
#' factor is the same as initial chain).
#' @param keep_old_samples Logical indicating if the original samples should be
#' kept or only the new samples returned. Defaults to TRUE.
#' @param ... Accepts only the deprecated \code{R} argument (renamed to
#' \code{n_iter}; still works, with a warning, for this release, when
#' passed positionally in its original slot or when every other argument is
#' also named - naming \code{R} while leaving later arguments positional is
#' not supported, since there is no way to bind two names to the same
#' argument slot). Anything else is an "unused argument" error.
#' @return A named list containing the sampled partitions, cluster and batch
#' parameters, model fit measures and some details on the model call.
#' @export
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
#' fixed <- c(rep(1, 10), rep(0, 40), rep(1, 10), rep(0, 40))
#'
#' # Batch
#' batch_vec <- sample(seq(1, 5), replace = TRUE, size = 100)
#'
#' # Density choice
#' type <- "MVT"
#'
#' # Sampling parameters
#' n_iter <- 1000
#' thin <- 50
#' n_chains <- 4
#'
#' # MCMC samples
#' mcmc_output <- fitBatchMix(
#'   X,
#'   n_chains,
#'   n_iter,
#'   thin,
#'   batch_vec,
#'   type,
#'   initial_labels = labels,
#'   fixed = fixed
#' )
#'
#' # Given an initial value for the parameters
#' new_output <- continueChains(
#'   mcmc_output,
#'   X,
#'   fixed,
#'   batch_vec,
#'   n_iter,
#'   keep_old_samples = TRUE
#' )
continueChains <- function(mcmc_output,
                           X,
                           fixed,
                           batch_vec,
                           n_iter,
                           keep_old_samples = TRUE,
                           ...) {
  n_iter <- .resolveDeprecatedNIter(
    missing(n_iter), if (missing(n_iter)) NULL else n_iter, list(...), "continueChains"
  )

  new_output <- lapply(
    mcmc_output,
    continueChain,
    X,
    fixed,
    batch_vec,
    n_iter,
    keep_old_samples
  )

  n_chains <- length(mcmc_output)

  # Record chain number
  for (ii in seq(n_chains)) {
    new_output[[ii]]$Chain <- mcmc_output[[ii]]$Chain
  }

  # Convergence (see ``runMCMCChains``) is recomputed, not copied from
  # ``mcmc_output``: the chains are now longer, so the old Rhat/ESS/best
  # chain are stale.
  if (n_chains >= 2) {
    convergence <- tryCatch(assessConvergence(new_output), error = function(e) NULL)
    if (!is.null(convergence)) {
      attr(new_output, "convergence") <- convergence
      attr(new_output, "best_chain") <- convergence$best_chain
    }
  }

  # lapply() (above) drops the batchmix_fit_list class of the input list;
  # each element is still a classed batchmix_fit though (continueChain()
  # builds its output from batchSemiSupervisedMixtureModel()), so restore it
  # on the outer list too - see R/batchmixFitMethods.R.
  class(new_output) <- c("batchmix_fit_list", class(new_output))

  new_output
}
