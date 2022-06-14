#' @title Continue chains
#' @description Continues sampling from a list of previous chains.
#' @param mcmc_output Chains to be continued.
#' @param X Data to cluster as a matrix with the items to cluster held in rows.
#' @param fixed The indicator vector for which labels are observed.
#' @param batch_vec The vector of the batch labels for the data.
#' @param R The number of iterations to run in this continuation (thinning
#' factor is the same as initial chain).
#' @param keep_old_samples Logical indicating if the original samples should be
#' kept or only the new samples returned. Defaults to TRUE.
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
#' R <- 1000
#' thin <- 50
#' n_chains <- 4
#'
#' # MCMC samples and BIC vector
#' mcmc_output <- runMCMCChains(
#'   X,
#'   n_chains,
#'   R,
#'   thin,
#'   labels,
#'   fixed,
#'   batch_vec,
#'   type
#' )
#'
#' # Given an initial value for the parameters
#' new_output <- continueChains(
#'   mcmc_output,
#'   X,
#'   fixed,
#'   batch_vec,
#'   R,
#'   keep_old_samples = TRUE
#' )
continueChains <- function(mcmc_output,
                           X,
                           fixed,
                           batch_vec,
                           R,
                           keep_old_samples = TRUE) {
  new_output <- lapply(
    mcmc_output,
    continueChain,
    X,
    fixed,
    batch_vec,
    R,
    keep_old_samples
  )
  
  n_chains <- length(mcmc_output)
  
  # Record chain number 
  for(ii in seq(n_chains)) {
    new_output[[ii]]$Chain <- mcmc_output[[ii]]$Chain
  }
  

  new_output
}
