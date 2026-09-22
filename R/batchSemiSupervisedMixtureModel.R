#' @title Batch semi-supervised/unsupervised mixture model (low-level engine)
#' @description The low-level engine that \code{\link{runBatchMix}} and
#' (through it) \code{\link{fitBatchMix}} build on - despite the name, it
#' handles both the semi-supervised (\code{fixed} has any 1s) and fully
#' unsupervised (\code{fixed} all 0, the default) case identically to those
#' wrappers. Most users should call \code{\link{fitBatchMix}} (fits multiple
#' chains, the recommended entry point) or \code{\link{runBatchMix}} (a
#' single chain, with convenience defaults for \code{K_max}/\code{initial_labels}
#' this function does not provide) instead of this function directly - see
#' \code{?runBatchMix} for the full, canonical description of every argument
#' below (inherited via \code{@@inheritParams}, not repeated here, so there is
#' exactly one place to read about them).
#' @inheritParams runBatchMix
#' @param concentration Initial concentration vector for component weights;
#' alternative to \code{alpha} (only one of the two should be given) for
#' directly specifying an asymmetric concentration per component.
#' @return A named list containing the sampled partitions, cluster and batch
#' parameters, model fit measures and some details on the model call.
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
#' type <- "MVN"
#'
#' # Sampling parameters
#' n_iter <- 1000
#' thin <- 50
#'
#' # MCMC samples and BIC vector
#' samples <- batchSemiSupervisedMixtureModel(
#'   X,
#'   n_iter,
#'   thin,
#'   labels,
#'   fixed,
#'   batch_vec,
#'   type
#' )
#'
#' # Given an initial value for the parameters
#' initial_class_means <- matrix(c(1, 1, 3, 4), nrow = 2)
#' initial_class_covariance <- array(c(1, 0, 0, 1, 1, 0, 0, 1),
#'   dim = c(2, 2, 2)
#' )
#'
#' # We can use values from a previous chain
#' initial_batch_shift <- samples$batch_shift[, , n_iter / thin]
#' initial_batch_scale <- matrix(
#'   c(1.2, 1.3, 1.7, 1.1, 1.4, 1.3, 1.2, 1.2, 1.1, 2.0),
#'   nrow = 2
#' )
#'
#' samples <- batchSemiSupervisedMixtureModel(X,
#'   n_iter,
#'   thin,
#'   labels,
#'   fixed,
#'   batch_vec,
#'   type,
#'   initial_class_means = initial_class_means,
#'   initial_class_covariance = initial_class_covariance,
#'   initial_batch_shift = initial_batch_shift,
#'   initial_batch_scale = initial_batch_scale
#' )
#'
batchSemiSupervisedMixtureModel <- function(X,
                                            n_iter,
                                            thin,
                                            initial_labels,
                                            fixed,
                                            batch_vec,
                                            type,
                                            # -- problem specification --
                                            K_max = length(unique(initial_labels)),
                                            alpha = NULL,
                                            concentration = NULL,
                                            # -- MCMC control (auto-tuning is on by default - see Description) --
                                            auto_tune = TRUE,
                                            n_burn = NULL,
                                            verbose = TRUE,
                                            # -- proposal windows (only matter if auto_tune = FALSE) --
                                            mu_proposal_window = 0.5**2,
                                            cov_proposal_window = 0.002,
                                            r_proposal_window = 0.1,
                                            sigma_proposal_window = 0.01,
                                            m_proposal_window = 0.3**2,
                                            S_proposal_window = 0.01,
                                            t_df_proposal_window = 0.015,
                                            gamma_proposal_window = 0.1,
                                            eta_proposal_window = 0.1,
                                            gp_hyperparameter_proposal_window = 0.1,
                                            # -- prior hyperparameters --
                                            m_scale = NULL,
                                            rho = 3.0,
                                            theta = 1.0,
                                            eta = 1.0,
                                            a_gamma = 2.0,
                                            b_gamma = 1.0,
                                            # -- initial values (warm starts; default to prior draws) --
                                            initial_class_means = NULL,
                                            initial_class_covariance = NULL,
                                            initial_batch_shift = NULL,
                                            initial_batch_scale = NULL,
                                            initial_class_df = NULL,
                                            # -- optional structural extensions --
                                            include_interaction = FALSE,
                                            batch_weight_prior = c("global", "partial_pooling", "gp"),
                                            batch_coordinates = NULL,
                                            gp_tau2 = 1.0,
                                            gp_length_scale = 1.0,
                                            sample_gp_hyperparameters = FALSE,
                                            pp_tau2_shape = 2.0,
                                            pp_tau2_rate = 1.0,
                                            pp_mu_prior_sd = 10.0,
                                            # -- MVN_MIXED-only --
                                            column_type = NULL,
                                            censor_code = NULL,
                                            ...) {
  n_iter <- .resolveDeprecatedNIter(
    missing(n_iter), if (missing(n_iter)) NULL else n_iter, list(...),
    "batchSemiSupervisedMixtureModel"
  )

  if (!is.matrix(X)) {
    stop("X is not a matrix. Data should be in matrix format.")
  }

  if (length(batch_vec) != nrow(X)) {
    stop("The number of rows in X and the number of batch labels are not equal.")
  }

  # Only the "MVN_MIXED" sampler (mvnSamplerMixed) implements missing-data
  # augmentation. The other three types set up their priors from raw
  # mean(X)/cov(X) with no NA-handling, so a missing entry silently poisons
  # every downstream calculation and the sampler will fail deep inside the
  # C++ layer (e.g. an "not symmetric positive definite" error from the
  # Inverse-Wishart prior) rather than at this, more informative, entry
  # point. Catch it here instead.
  if (type != "MVN_MIXED" && anyNA(X)) {
    stop(paste0(
      "X contains missing values (NA/NaN), but type = '", type, "' does not ",
      "support missing data. Use type = 'MVN_MIXED' instead, which models ",
      "missing (and censored) entries via proper data augmentation - see ",
      "?sampleSemisupervisedMVNMixed and the 'probit_missing_censored' vignette."
    ))
  }

  if (n_iter < thin) {
    warning("Iterations to run less than thinning factor. No samples recorded.")
  }

  # Check that the initial labels starts at 0, if not remedy this.
  if (!any(initial_labels == 0)) {
    initial_labels <- as.numeric(as.factor(initial_labels)) - 1
  }

  if (max(initial_labels) != (length(unique(initial_labels)) - 1)) {
    stop("initial labels are not all contiguous integers.")
  }

  # Check that the batch labels starts at 0, if not remedy this.
  if (!any(batch_vec == 0)) {
    batch_vec <- as.numeric(as.factor(batch_vec)) - 1
  }

  if (max(batch_vec) != (length(unique(batch_vec)) - 1)) {
    stop("batch labels are not all contiguous integers.")
  }

  # The number of batches present
  B <- length(unique(batch_vec))

  # The concentration parameter for the prior Dirichlet distribution of the
  # component weights.
  alpha_not_passed <- is.null(alpha)
  concentration_not_passed <- is.null(concentration)
  if ((!alpha_not_passed) && (!concentration_not_passed)) {
    stop("Only one of ``concentration`` or ``alpha`` should be passed.")
  }
  if (concentration_not_passed) {
    if (alpha_not_passed) {
      alpha <- 1.0 / K_max
    }
    concentration <- rep(alpha, K_max)
  }

  # Check the proposal windows are all strictly positive
  checkProposalWindows(
    mu_proposal_window,
    cov_proposal_window,
    m_proposal_window,
    S_proposal_window,
    t_df_proposal_window,
    verbose
  )

  # The proposal windows for these objects are narrower for larger quantities,
  # so we use the reciprocal to ensure that the relationship between acceptance
  # rates is the same for all parameters, namely that smaller windows increases
  # the acceptance rate
  actual_cov_proposal_window <- 1.0 / cov_proposal_window
  actual_S_proposal_window <- 1.0 / S_proposal_window
  actual_t_df_proposal_window <- 1.0 / t_df_proposal_window
  actual_sigma_proposal_window <- 1.0 / sigma_proposal_window

  P <- ncol(X)

  # Auto-tuning defaults to adapting over the first half of the run.
  if (is.null(n_burn)) {
    n_burn <- floor(n_iter / 2)
  }

  # Validate and map to the integer code the C++ layer expects: 0 = global,
  # 1 = partial pooling, 2 = gp.
  batch_weight_prior <- match.arg(batch_weight_prior)
  weight_prior_type <- switch(batch_weight_prior,
    global = 0L,
    partial_pooling = 1L,
    gp = 2L
  )

  # The 1-D coordinate for each batch, used only if batch_weight_prior is
  # "gp"; defaults to batch collection/label order.
  if (is.null(batch_coordinates)) {
    batch_coordinates <- numeric(0)
  } else if (length(batch_coordinates) != B) {
    stop("batch_coordinates must have one entry per batch (length B).")
  }

  # Only meaningful for type = 'MVN_MIXED'.
  if (is.null(column_type)) {
    column_type <- rep(0L, P)
  } else if (length(column_type) != P) {
    stop("column_type must have one entry per column of X (length P).")
  }
  if (is.null(censor_code)) {
    censor_code <- matrix(0L, nrow = nrow(X), ncol = P)
  } else if (!all(dim(censor_code) == c(nrow(X), P))) {
    stop("censor_code must be an N x P matrix matching the dimensions of X.")
  }

  # Check if an initial value is passed for any of the parameters. Prepare the
  # parameters to be passed to C++.
  class_mean_passed <- !is.null(initial_class_means)
  class_covariance_passed <- !is.null(initial_class_covariance)
  batch_shift_passed <- !is.null(initial_batch_shift)
  batch_scale_passed <- !is.null(initial_batch_scale)
  class_df_passed <- !is.null(initial_class_df)

  sample_m_scale <- is.null(m_scale)
  if(sample_m_scale) {
    m_scale <- 0.01
  }

  initial_parameters <- prepareInitialParameters(
    initial_class_means,
    initial_class_covariance,
    initial_batch_shift,
    initial_batch_scale,
    initial_class_df,
    P,
    K_max,
    B,
    type
  )

  class_means <- initial_parameters$class_means
  class_cov <- initial_parameters$class_cov
  batch_shift <- initial_parameters$batch_shift
  batch_scale <- initial_parameters$batch_scale
  class_df <- initial_parameters$class_df

  # Common auto-tuning / interaction / GP-correlated-weight arguments shared
  # by every sampler type.
  extra_args <- list(
    auto_tune = auto_tune,
    n_burn = n_burn,
    include_interaction = include_interaction,
    gamma_proposal_window = gamma_proposal_window,
    a_gamma = a_gamma,
    b_gamma = b_gamma,
    weight_prior_type = weight_prior_type,
    batch_coordinates = batch_coordinates,
    gp_tau2 = gp_tau2,
    gp_length_scale = gp_length_scale,
    eta_proposal_window = eta_proposal_window,
    sample_gp_hyperparameters = sample_gp_hyperparameters,
    gp_hyperparameter_proposal_window = gp_hyperparameter_proposal_window,
    pp_tau2_shape = pp_tau2_shape,
    pp_tau2_rate = pp_tau2_rate,
    pp_mu_prior_sd = pp_mu_prior_sd
  )

  # Pull samples from the mixture model
  if (type == "MVN") {
    mcmc_output <- do.call(sampleSemisupervisedMVN, c(list(
      X,
      K_max,
      B,
      initial_labels,
      batch_vec,
      fixed,
      mu_proposal_window,
      actual_cov_proposal_window,
      m_proposal_window,
      actual_S_proposal_window,
      n_iter,
      thin,
      concentration,
      m_scale,
      rho,
      theta,
      class_means,
      class_cov,
      batch_shift,
      batch_scale,
      class_mean_passed,
      class_covariance_passed,
      TRUE, # batch_shift_passed,
      batch_scale_passed,
      sample_m_scale
    ), extra_args))
  }

  if (type == "MVT") {
    mcmc_output <- do.call(sampleSemisupervisedMVT, c(list(
      X,
      K_max,
      B,
      initial_labels,
      batch_vec,
      fixed,
      mu_proposal_window,
      actual_cov_proposal_window,
      m_proposal_window,
      actual_S_proposal_window,
      actual_t_df_proposal_window,
      n_iter,
      thin,
      concentration,
      m_scale,
      rho,
      theta,
      class_means,
      class_cov,
      class_df,
      batch_shift,
      batch_scale,
      class_mean_passed,
      class_covariance_passed,
      class_df_passed,
      TRUE, # batch_shift_passed,
      batch_scale_passed,
      sample_m_scale
    ), extra_args))
  }

  if (type == "MVN_LKJ") {
    mcmc_output <- do.call(sampleSemisupervisedMVNSeparationStrategy, c(list(
      X,
      K_max,
      B,
      initial_labels,
      batch_vec,
      fixed,
      mu_proposal_window,
      r_proposal_window,
      actual_sigma_proposal_window,
      m_proposal_window,
      actual_S_proposal_window,
      n_iter,
      thin,
      concentration,
      m_scale,
      rho,
      theta,
      class_means,
      class_cov,
      batch_shift,
      batch_scale,
      class_mean_passed,
      class_covariance_passed,
      TRUE, # batch_shift_passed,
      batch_scale_passed,
      sample_m_scale,
      eta
    ), extra_args))
  }

  if (type == "MVN_MIXED") {
    mcmc_output <- do.call(sampleSemisupervisedMVNMixed, c(list(
      X,
      K_max,
      B,
      initial_labels,
      batch_vec,
      fixed,
      column_type,
      censor_code,
      mu_proposal_window,
      r_proposal_window,
      actual_sigma_proposal_window,
      m_proposal_window,
      actual_S_proposal_window,
      n_iter,
      thin,
      concentration,
      m_scale,
      rho,
      theta,
      eta,
      sample_m_scale
    ), extra_args))
  }

  mcmc_output$concentration <- matrix(concentration, nrow = 1)


  if (!type %in% c("MVN", "MVT", "MVN_LKJ", "MVN_MIXED")) {
    stop("Type not recognised. Please use one of 'MVN', 'MVT', 'MVN_LKJ' or 'MVN_MIXED'.")
  }

  # Record details of model run to output
  # MCMC details
  mcmc_output$thin <- thin
  mcmc_output$n_iter <- n_iter
  mcmc_output$burn <- 0

  # Density choice
  mcmc_output$type <- type

  # Dimensions of data
  mcmc_output$P <- P
  mcmc_output$N <- nrow(X)

  # Number of components and batches modelled
  mcmc_output$K_max <- K_max
  mcmc_output$B <- B

  # Record hyperparameter choice
  mcmc_output$alpha <- alpha
  mcmc_output$m_scale <- m_scale
  mcmc_output$rho <- rho
  mcmc_output$theta <- theta

  # Proposal windows
  mcmc_output$mu_proposal_window <- mu_proposal_window
  mcmc_output$cov_proposal_window <- cov_proposal_window
  mcmc_output$r_proposal_window <- r_proposal_window
  mcmc_output$sigma_proposal_window <- sigma_proposal_window
  mcmc_output$m_proposal_window <- m_proposal_window
  mcmc_output$S_proposal_window <- S_proposal_window
  mcmc_output$t_df_proposal_window <- t_df_proposal_window

  # Auto-tuning / interaction / GP-correlated-weight model choices
  mcmc_output$auto_tune <- auto_tune
  mcmc_output$n_burn <- n_burn
  mcmc_output$include_interaction <- include_interaction
  mcmc_output$gamma_proposal_window <- gamma_proposal_window
  mcmc_output$a_gamma <- a_gamma
  mcmc_output$b_gamma <- b_gamma
  mcmc_output$batch_weight_prior <- batch_weight_prior
  mcmc_output$batch_coordinates <- batch_coordinates
  mcmc_output$eta_proposal_window <- eta_proposal_window
  mcmc_output$sample_gp_hyperparameters <- sample_gp_hyperparameters
  mcmc_output$gp_hyperparameter_proposal_window <- gp_hyperparameter_proposal_window
  mcmc_output$pp_tau2_shape <- pp_tau2_shape
  mcmc_output$pp_tau2_rate <- pp_tau2_rate
  mcmc_output$pp_mu_prior_sd <- pp_mu_prior_sd
  mcmc_output$eta <- eta

  # Indicate if the model was semi-supervised or unsupervised
  mcmc_output$Semisupervised <- TRUE

  # Correct this if we were effectively unsupervised
  actually_unsupervised <- all(fixed == 0)
  if (actually_unsupervised) {
    mcmc_output$Semisupervised <- FALSE
  }
  
  # Indicate if lambda^2 was sampled
  mcmc_output$sample_m_scale <- sample_m_scale
  
  mcmc_output
}
