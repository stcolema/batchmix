#!/usr/bin/Rscript
#' @title Collect acceptance rate
#' @description Collects the acceptance rates for each parameter into a data.frame
#' @param samples The output of ``runBatchMix``.
#' @return A wide data.frame of all the sampled parameters and the iteration.
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
#'
#' # MCMC samples
#' samples <- runBatchMix(X, R, thin, batch_vec, "MVN",
#'   initial_labels = labels,
#'   fixed = fixed
#' )
#'
#' # Acceptance rates
#' collectAcceptanceRates(samples)
#'
collectAcceptanceRates <- function(samples) {
  # Number of classes and batches
  K <- samples$K_max # ncol(samples$means[, , 1])
  B <- samples$B # ncol(samples$batch_shift[, , 1])
  P <- samples$P # nrow(samples$means[, , 1])
  N <- samples$N # ncol(samples$samples)
  type <- samples$type

  # Stack the sampled matrices on top of each other
  mean_rate <- t(samples$mu_acceptance_rate)
  colnames(mean_rate) <- paste0("Mu_", seq(1, K))

  # The covariance parameterisation (and hence its acceptance-rate fields)
  # differs by type: MVN/MVT use a single Wishart-proposed covariance,
  # MVN_LKJ/MVN_MIXED use the separate R (correlation) and sigma (marginal
  # scale) proposals.
  if (type %in% c("MVN", "MVT")) {
    cov_rate <- t(samples$cov_acceptance_rate)
    colnames(cov_rate) <- paste0("Sigma_", seq(1, K))
  } else {
    r_rate <- t(samples$r_acceptance_rate)
    colnames(r_rate) <- paste0("R_", seq(1, K))
    sigma_rate <- t(samples$sigma_acceptance_rate)
    colnames(sigma_rate) <- paste0("sigma_", seq(1, K))
    cov_rate <- cbind(r_rate, sigma_rate)
  }

  batch_shift_rate <- t(samples$m_acceptance_rate)
  colnames(batch_shift_rate) <- paste0("m_", seq(1, B))
  batch_scale_rate <- t(samples$S_acceptance_rate)
  colnames(batch_scale_rate) <- paste0("S_", seq(1, B))

  output_df <- as.data.frame(
    cbind(
      mean_rate,
      cov_rate,
      batch_shift_rate,
      batch_scale_rate
    )
  )

  if (type == "MVT") {
    t_df_rate <- t(samples$t_df_acceptance_rate)
    colnames(t_df_rate) <- paste0("nu_", seq(1, K))
    output_df <- cbind(output_df, t_df_rate)
  }

  if (type == "MSN") {
    phi_rate <- t(samples$phi_acceptance_rate)
    colnames(phi_rate) <- paste0("phi_", seq(1, K))
    output_df <- cbind(output_df, phi_rate)
  }

  # Interaction term acceptance rate (one entry per feature dimension P;
  # see batchSemiSupervisedMixtureModel()'s include_interaction argument).
  if (isTRUE(samples$include_interaction) && !is.null(samples$gamma_acceptance_rate)) {
    gamma_rate <- matrix(samples$gamma_acceptance_rate, nrow = 1)
    colnames(gamma_rate) <- paste0("gamma_", seq_len(ncol(gamma_rate)))
    output_df <- cbind(output_df, gamma_rate)
  }

  # Batch-specific weight acceptance rates (one entry per free
  # additive-log-ratio coordinate, K - 1 of them; see the
  # batch_weight_prior argument - "partial_pooling" or "gp").
  batch_weights_used <- isTRUE(samples$weight_prior_type > 0) || isTRUE(samples$batch_weight_prior %in% c("partial_pooling", "gp"))
  if (batch_weights_used && !is.null(samples$eta_acceptance_rate)) {
    eta_rate <- matrix(samples$eta_acceptance_rate, nrow = 1)
    colnames(eta_rate) <- paste0("eta_", seq_len(ncol(eta_rate)))
    output_df <- cbind(output_df, eta_rate)
    if (!is.null(samples$gp_hyperparameter_acceptance_rate)) {
      output_df$gp_hyperparameter <- samples$gp_hyperparameter_acceptance_rate
    }
  }

  output_df
}
