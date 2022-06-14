#!/usr/bin/Rscript
#' @title Run MCMC Chains
#' @description Run multiple chains of the semi-supervised batch mixture model
#' of the same type and from the same initial clustering.
#' @param X Data to cluster as a matrix with the items to cluster held in rows.
#' @param n_chains Integer. Number of MCMC chains to run.
#' @param initial_labels Initial clustering.
#' @param fixed Which items are fixed in their initial label.
#' @param batch_vec Labels identifying which batch each item being clustered is 
#' from.
#' @param R The number of iterations in the sampler.
#' @param thin The factor by which the samples generated are thinned, e.g. if
#' ``thin=50`` only every 50th sample is kept.
#' @param type Character indicating density type to use. One of 'MVN'
#' (multivariate normal distribution), 'MVT' (multivariate t distribution) or
#' 'MSN' (multivariate skew normal distribution).
#' @param K_max The number of components to include (the upper bound on the
#' number of clusters in each sample). Defaults to the number of unique labels
#' in ``initial_labels``.
#' @param alpha The concentration parameter for the stick-breaking prior and the
#' weights in the model.
#' @param mu_proposal_window The proposal window for the cluster mean proposal
#' kernel.
#' @param cov_proposal_window The proposal window for the cluster covariance
#' proposal kernel.
#' @param m_proposal_window The proposal window for the batch mean proposal
#'  kernel.
#' @param S_proposal_window The proposal window for the batch standard deviation
#'  proposal kernel.
#' @param t_df_proposal_window The proposal window for the degrees of freedom
#' for the multivariate t distribution (not used if type is not 'MVT').
#' @param phi_proposal_window The proposal window for the shape parameter for
#' the multivariate skew normal distribution (not used if type is not 'MSN').
#' @param m_scale The scale hyperparameter for the batch shift prior
#' distribution.
#' @param rho The shape of the prior distribution for the batch scale.
#' @param theta The scale of the prior distribution for the batch scale.
#' @returns A list of named lists. Each entry is the output of 
#' ``batchSemiSupervisedMixtureModel``.
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
#' # MCMC samples and BIC vector
#' samples <- runMCMCChains(X, n_chains, R, thin, labels, fixed, batch_vec, "MVN")
#' 
runMCMCChains <- function(X,
                          n_chains,
                          R,
                          thin,
                          initial_labels,
                          fixed,
                          batch_vec,
                          type,
                          K_max = length(unique(initial_labels)),
                          alpha = NULL,
                          mu_proposal_window = 0.5**2,
                          cov_proposal_window = 100,
                          m_proposal_window = 0.3**2,
                          S_proposal_window = 100,
                          t_df_proposal_window = 100,
                          phi_proposal_window = 1.2**2,
                          m_scale = 0.01,
                          rho = 3.0,
                          theta = 1.0) {
  mcmc_lst <- vector("list", n_chains)

  mcmc_lst <- lapply(mcmc_lst, function(x) {
    batchSemiSupervisedMixtureModel(X,
      R,
      thin,
      initial_labels,
      fixed,
      batch_vec,
      type,
      K_max = K_max,
      alpha = alpha,
      mu_proposal_window = mu_proposal_window,
      cov_proposal_window = cov_proposal_window,
      m_proposal_window = m_proposal_window,
      S_proposal_window = S_proposal_window,
      t_df_proposal_window = t_df_proposal_window,
      phi_proposal_window = phi_proposal_window,
      m_scale = m_scale,
      rho = rho,
      theta = theta
    )
  })
  
  # Record chain number 
  for(ii in seq(n_chains)) {
    mcmc_lst[[ii]]$Chain <- ii
  }
  
  mcmc_lst
}
