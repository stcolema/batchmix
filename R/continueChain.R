#' @title Continue chain
#' @description Continues sampling from a previous position for a given chain.
#' @param mcmc_output Chain to be continued.
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
#'
#' # MCMC samples and BIC vector
#' mcmc_output <- batchSemiSupervisedMixtureModel(
#'   X,
#'   R,
#'   thin,
#'   labels,
#'   fixed,
#'   batch_vec,
#'   type
#' )
#'
#' # Given an initial value for the parameters
#' mcmc_output <- continueChain(
#'   mcmc_output,
#'   X,
#'   fixed,
#'   batch_vec,
#'   R,
#' )
#' 
continueChain <- function(mcmc_output,
                          X,
                          fixed,
                          batch_vec,
                          R,
                          keep_old_samples = TRUE) {

  # The relevant aspects of the previous chain
  R_old <- mcmc_output$R
  thin <- mcmc_output$thin
  last_sample <- R_eff_old <- floor(R_old / thin)

  B <- mcmc_output$B
  K_max <- mcmc_output$K_max
  N <- mcmc_output$N
  P <- mcmc_output$P

  type <- mcmc_output$type

  alpha <- mcmc_output$alpha
  m_scale <- mcmc_output$m_scale
  rho <- mcmc_output$rho
  theta <- mcmc_output$theta

  mu_proposal_window <- mcmc_output$mu_proposal_window
  cov_proposal_window <- mcmc_output$cov_proposal_window
  m_proposal_window <- mcmc_output$m_proposal_window
  S_proposal_window <- mcmc_output$S_proposal_window
  t_df_proposal_window <- mcmc_output$t_df_proposal_window
  phi_proposal_window <- mcmc_output$phi_proposal_window

  initial_class_means <- mcmc_output$means[, , last_sample]

  initial_class_covariance <- array(0, dim = c(P, P, K_max))
  for (k in seq(1, K_max)) {
    col_i <- (k - 1) * P + 1
    col_j <- k * P
    rel_cols <- seq(col_i, col_j)
    initial_class_covariance[, , k] <- mcmc_output$covariance[, rel_cols, last_sample]
  }

  initial_batch_shift <- mcmc_output$batch_shift[, , last_sample]
  initial_batch_scale <- mcmc_output$batch_scale[, , last_sample]

  is_mvt <- type == "MVT"
  is_semisupervised <- mcmc_output$Semisupervised

  if (is_mvt) {
    initial_class_df <- mcmc_output$t_df[last_sample, ]
  }


  labels <- mcmc_output$samples[last_sample, ]

  new_samples <- batchSemiSupervisedMixtureModel(X,
    R,
    thin,
    labels,
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
    theta = theta,
    initial_class_means = initial_class_means,
    initial_class_covariance = initial_class_covariance,
    initial_class_df = initial_class_df,
    initial_batch_shift = initial_batch_shift,
    initial_batch_scale = initial_batch_scale
  )

  if (keep_old_samples) {
    R_eff <- floor(R / thin)

    R_comb_eff <- R_eff + R_eff_old

    # R_comb is has to be calculated this way. Reason, consider R_old = 107,
    # R = 113, thin = 10; then R_old + R = 220, 220 / 10 = 22, but we only
    # record R_eff_old = floor(107 / 10) = 10, R_eff = floor(113 / 10) = 11,
    # R_comb_eff = 21.
    R_comb <- R_comb_eff * thin

    # Let's combine the sampled parameters
    combined_allocation_samples <- rbind(
      mcmc_output$samples,
      new_samples$samples
    )

    combined_means <- array(c(mcmc_output$means, new_samples$means),
      dim = c(P, K_max, R_comb_eff)
    )

    combined_shifts <- array(
      c(
        mcmc_output$batch_shift,
        new_samples$batch_shift
      ),
      dim = c(P, B, R_comb_eff)
    )

    combined_scales <- array(
      c(
        mcmc_output$batch_scale,
        new_samples$batch_scale
      ),
      dim = c(P, B, R_comb_eff)
    )

    combined_mean_sums <- array(
      c(
        mcmc_output$mean_sum,
        new_samples$mean_sum
      ),
      dim = c(P, K_max * B, R_comb_eff)
    )

    combined_covariances <- array(c(
      mcmc_output$covariance,
      new_samples$covariance
    ),
    dim = c(P, P * K_max, R_comb_eff)
    )

    combined_covariance_comb <- array(
      c(
        mcmc_output$cov_comb,
        new_samples$cov_comb
      ),
      dim = c(P, P * K_max * B, R_comb_eff)
    )

    combined_weights <- rbind(mcmc_output$weights, new_samples$weights)

    comb_cov_acceptance_rate <- ((mcmc_output$cov_acceptance_rate * R_old +
      new_samples$cov_acceptance_rate * R)
    / (R_old + R)
    )

    comb_mu_acceptance_rate <- ((mcmc_output$mu_acceptance_rate * R_old +
      new_samples$mu_acceptance_rate * R)
    / (R_old + R)
    )

    comb_m_acceptance_rate <- ((mcmc_output$m_acceptance_rate * R_old +
      new_samples$m_acceptance_rate * R)
    / (R_old + R)
    )

    comb_S_acceptance_rate <- ((mcmc_output$S_acceptance_rate * R_old +
      new_samples$S_acceptance_rate * R)
    / (R_old + R)
    )

    if (type == "MVT") {
      comb_t_df_acceptance_rate <- ((mcmc_output$t_df_acceptance_rate * R_old +
        new_samples$t_df_acceptance_rate * R)
      / (R_old + R)
      )

      comb_t_df <- rbind(mcmc_output$t_df, new_samples$t_df)
    }

    # if(is_semisupervised) {
    comb_allocation_probs <- array(c(
      mcmc_output$alloc,
      new_samples$alloc
    ),
    dim = c(N, K_max, R_comb_eff)
    )

    # }

    comb_complete_likelihood <- rbind(
      mcmc_output$complete_likelihood,
      new_samples$complete_likelihood
    )

    comb_observed_likelihood <- rbind(
      mcmc_output$observed_likelihood,
      new_samples$observed_likelihood
    )

    comb_bic <- rbind(
      mcmc_output$BIC,
      new_samples$BIC
    )


    comb_inferred_data <- array(c(
      mcmc_output$batch_corrected_data,
      new_samples$batch_corrected_data
    ),
    dim = c(N, P, R_comb_eff)
    )

    R_comb

    new_samples$R <- R_comb

    new_samples$means <- combined_means
    new_samples$covariance <- combined_covariances
    new_samples$batch_shift <- combined_shifts
    new_samples$batch_scale <- combined_scales
    new_samples$mean_sum <- combined_mean_sums
    new_samples$cov_comb <- combined_covariance_comb

    new_samples$weights <- combined_weights
    new_samples$cov_acceptance_rate <- comb_cov_acceptance_rate
    new_samples$mu_acceptance_rate <- comb_mu_acceptance_rate
    new_samples$m_acceptance_rate <- comb_m_acceptance_rate
    new_samples$S_acceptance_rate <- comb_S_acceptance_rate

    new_samples$observed_likelihood <- comb_observed_likelihood
    new_samples$complete_likelihood <- comb_complete_likelihood
    new_samples$BIC <- comb_bic

    new_samples$samples <- combined_allocation_samples
    new_samples$alloc <- comb_allocation_probs
    new_samples$batch_corrected_data <- comb_inferred_data

    if (is_mvt) {
      new_samples$t_df <- comb_t_df
      new_samples$t_df_acceptance_rate <- comb_t_df_acceptance_rate
    }
  }

  new_samples
}
