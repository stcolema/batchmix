#!/usr/bin/Rscript
#' @title Run MCMC Chains
#' @description Run multiple chains of the batch mixture model of the same type.
#' @param X Data to cluster as a matrix with the items to cluster held in rows.
#' @param n_chains Integer. Number of MCMC chains to run.
#' @param batch_vec Labels identifying which batch each item being clustered is
#' from.
#' @param R The number of iterations in the sampler.
#' @param thin The factor by which the samples generated are thinned, e.g. if
#' ``thin=50`` only every 50th sample is kept.
#' @param type Character indicating density type to use. One of 'MVN'
#' (multivariate normal distribution) or 'MVT' (multivariate t distribution).
#' @param batch_specific_weights Allow each batch to have unique class weights.
#' If FALSE class weights are common across batches, if TRUE the class weights
#' are nested within batch and share a common hyperparameter. Defaults to FALSE.
#' If batch-specific weights are used then the returned object contains a member
#' ``weights`` which is a matrix with K x B columns. The columns are ordered by
#' batch, i.e. the first K columns contain the class weights in the first batch,
#' the second K are the class weights in the second batch, etc. If generic
#' weights are used then this matrix has K columns, one for each component weight.
#' @param K_max The number of components to include (the upper bound on the
#' number of clusters in each sample). Defaults to the number of unique labels
#' in ``initial_labels``.
#' @param initial_labels Initial clustering, if none given defaults to a random draw.
#' @param fixed Which items are fixed in their initial label. If not given,
#' defaults to a vector of 0 meaning the model is run unsupervised.
#' @param alpha The concentration parameter for the stick-breaking prior and the
#' weights in the model.
#' @param mu_proposal_window The proposal window for the cluster mean proposal
#' kernel. The proposal density is a Gaussian distribution, the window is the
#' variance.
#' @param cov_proposal_window The proposal window for the cluster covariance
#' proposal kernel. The proposal density is a Wishart distribution, this
#' argument is the reciprocal of the degree of freedom.
#' @param m_proposal_window The proposal window for the batch mean proposal
#'  kernel. The proposal density is a Gaussian distribution, the window is the
#' variance.
#' @param S_proposal_window The proposal window for the batch standard deviation
#'  proposal kernel. The proposal density is a Gamma distribution, this
#' argument is the reciprocal of the rate.
#' @param t_df_proposal_window The proposal window for the degrees of freedom
#' for the multivariate t distribution (not used if type is not 'MVT'). The
#' proposal density is a Gamma distribution, this argument is the reciprocal of
#' the rate.
#' @param m_scale The scale hyperparameter for the batch shift prior
#' distribution. This defines the scale of the batch effect upon the mean and
#' should be in (0, 1].
#' @param rho The shape of the prior distribution for the batch scale.
#' @param theta The scale of the prior distribution for the batch scale.
#' @param initial_class_means A $P x K$ matrix of initial values for the class
#' means. Defaults to draws from the prior distribution.
#' @param initial_class_covariance A $P x P x K$ array of initial values for
#' the class covariance matrices. Defaults to draws from the prior distribution.
#' @param initial_batch_shift A $P x B$ matrix of initial values for the batch
#' shift effect Defaults to draws from the prior distribution.
#' @param initial_batch_scale A $P x B$ matrix of initial values for the batch
#' scales Defaults to draws from the prior distribution.
#' @param initial_class_df A $K$ vector of initial values for the class degrees
#' of freedom. Defaults to draws from the prior distribution.
#' @param verbose Logiccal indicating if warning about proposal windows should
#' be printed.
#' @returns A list of named lists. Each entry is the output of
#' ``runBatchMix``.
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
#' # Sampling parameters
#' R <- 1000
#' thin <- 50
#' n_chains <- 4
#'
#' # MCMC samples
#' samples <- runMCMCChains(X, n_chains, R, thin, batch_vec, "MVN",
#'   initial_labels = labels,
#'   fixed = fixed
#' )
#'
runMCMCChains <- function(X,
                          n_chains,
                          R,
                          thin,
                          batch_vec,
                          type,
                          batch_specific_weights = FALSE,
                          K_max = NULL,
                          initial_labels = NULL,
                          fixed = NULL,
                          alpha = 1,
                          mu_proposal_window = 0.5**2,
                          cov_proposal_window = 0.002,
                          m_proposal_window = 0.3**2,
                          S_proposal_window = 0.01,
                          t_df_proposal_window = 0.015,
                          m_scale = 0.01,
                          rho = 3.0,
                          theta = 1.0,
                          initial_class_means = NULL,
                          initial_class_covariance = NULL,
                          initial_batch_shift = NULL,
                          initial_batch_scale = NULL,
                          initial_class_df = NULL,
                          verbose = TRUE) {
  mcmc_lst <- vector("list", n_chains)

  mcmc_lst <- lapply(mcmc_lst, function(x) {
    runBatchMix(X,
      R,
      thin,
      batch_vec,
      type,
      K_max = K_max,
      initial_labels = initial_labels,
      fixed = fixed,
      alpha = alpha,
      batch_specific_weights = batch_specific_weights,
      mu_proposal_window = mu_proposal_window,
      cov_proposal_window = cov_proposal_window,
      m_proposal_window = m_proposal_window,
      S_proposal_window = S_proposal_window,
      t_df_proposal_window = t_df_proposal_window,
      m_scale = m_scale,
      rho = rho,
      theta = theta,
      initial_class_means = initial_class_means,
      initial_class_covariance = initial_class_covariance,
      initial_batch_shift = initial_batch_shift,
      initial_batch_scale = initial_batch_scale,
      initial_class_df = initial_class_df,
      verbose = verbose
    )
  })

  # Record chain number
  for (ii in seq(n_chains)) {
    mcmc_lst[[ii]]$Chain <- ii
  }

  mcmc_lst
}
