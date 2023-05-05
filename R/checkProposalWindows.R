#!/usr/bin/Rscript
#' @title Check proposal windows
#' @description Checks the proposal windows are acceptable.
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
#' @param verbose Logical indicating if a warning should be printed if proposal
#' windows are outside their expected scale.
#' @return No return value, called for side effects
#' @examples
#' checkProposalWindows(0.1, 0.2, 0.3, 0.1, 0.4, 0.3)
checkProposalWindows <- function(mu_proposal_window,
                                 cov_proposal_window,
                                 m_proposal_window,
                                 S_proposal_window,
                                 t_df_proposal_window,
                                 verbose = TRUE) {
  windows <- c(
    mu_proposal_window,
    m_proposal_window,
    cov_proposal_window,
    S_proposal_window,
    t_df_proposal_window
  )

  windows_reciprocals_used <- c(
    cov_proposal_window,
    S_proposal_window,
    t_df_proposal_window
  )

  all_windows_positive <- all(windows > 0)
  if (!all_windows_positive) {
    stop(paste0(
      "All proposal windowws must be positive numbers. Normally they",
      " take values between 0 and 1."
    ))
  }

  unexpected_value_taken <- any(windows_reciprocals_used >= 1.0)
  if (unexpected_value_taken && verbose) {
    warning(paste0(
      "The proposal windows for the covariance matrix, the batch scale, the ",
      "degrees of freedom and the skew are normally expected to be less than 1.0."
    ))
  }
  NULL
}
