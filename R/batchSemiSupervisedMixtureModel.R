#' @title Batch semisupervised mixture model
#' @description A Bayesian mixture model with batch effects.
#' @param X Data to cluster as a matrix with the items to cluster held in rows.
#' @param initial_labels Initial clustering.
#' @param fixed Which items are fixed in their initial label.
#' @param batch_vec Labels identifying which batch each item being clustered is
#' from.
#' @param R The number of iterations in the sampler.
#' @param thin The factor by which the samples generated are thinned, e.g. if
#' ``thin=50`` only every 50th sample is kept.
#' @param type Character indicating density type to use. One of 'MVN'
#' (multivariate normal distribution, Inverse-Wishart covariance prior),
#' 'MVT' (multivariate t distribution), 'MVN_LKJ' (multivariate normal with
#' an LKJ prior on the cluster correlation structure and log-normal marginal
#' scales, decoupling correlation and scale beliefs - see
#' \code{vignette("lkj_correlation_recovery", package = "batchmix")}), or
#' 'MVN_MIXED' (as 'MVN_LKJ', extended to support binary/probit columns and
#' missing/censored continuous entries - see
#' \code{vignette("probit_missing_censored", package = "batchmix")}, and the
#' \code{column_type}/\code{censor_code} arguments).
#' @param K_max The number of components to include (the upper bound on the
#' number of clusters in each sample). Defaults to the number of unique labels
#' in ``initial_labels``.
#' @param alpha The concentration parameter for the stick-breaking prior and the
#' weights in the model.
#' @param concentration Initial concentration vector for component weights.
#' @param mu_proposal_window The proposal window for the cluster mean proposal
#' kernel. The proposal density is a Gaussian distribution, the window is the
#' variance.
#' @param cov_proposal_window The proposal window for the cluster covariance
#' proposal kernel when \code{type} is 'MVN' or 'MVT'. The proposal density
#' is a Wishart distribution, this argument is the reciprocal of the degree
#' of freedom.
#' @param r_proposal_window Only used if \code{type} is 'MVN_LKJ' or
#' 'MVN_MIXED': the standard deviation of the (unconstrained-space) Gaussian
#' random walk proposal for the cluster correlation matrix R. Smaller values
#' give a tighter proposal (higher acceptance, smaller steps).
#' @param sigma_proposal_window Only used if \code{type} is 'MVN_LKJ' or
#' 'MVN_MIXED': the proposal window for the cluster marginal standard
#' deviations. As with \code{S_proposal_window}, the proposal density is a
#' Gamma distribution and this argument is the reciprocal of the rate.
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
#' @param eta The LKJ concentration parameter for the cluster correlation
#' prior; only used if \code{type} is 'MVN_LKJ' or 'MVN_MIXED'. eta = 1 is
#' uniform over the space of correlation matrices, eta > 1 shrinks
#' correlations towards 0.
#' @param column_type Only used if \code{type} is 'MVN_MIXED': a P-vector, 0
#' for a continuous column, 1 for a binary column observed via a probit
#' link. Defaults to all-continuous.
#' @param censor_code Only used if \code{type} is 'MVN_MIXED': an N x P
#' matrix, meaningful only for continuous columns: 0 = not censored,
#' 1 = left-censored (true value below the recorded X entry), 2 =
#' right-censored (true value above the recorded X entry). Defaults to no
#' censoring. Missing entries are indicated by \code{NA}/\code{NaN} in
#' \code{X} directly, for both continuous and binary columns.
#' @param auto_tune Logical; if \code{TRUE} (the default), every proposal
#' window is adapted during the first \code{n_burn} iterations via
#' Robbins-Monro diminishing adaptation, instead of staying fixed at the
#' value passed in. Adaptation is frozen after \code{n_burn} iterations so
#' the post-burn-in chain retains the correct stationary distribution.
#' @param n_burn Number of iterations treated as burn-in for proposal-window
#' adaptation; ignored if \code{auto_tune} is \code{FALSE}. Defaults to half
#' of \code{R}.
#' @param include_interaction Logical; if \code{TRUE}, add a batch x cluster
#' interaction term to the mean, \eqn{\gamma_{k,b} \sim N(0,
#' \tau^2_{\mathrm{interaction}})}, with \eqn{\tau^2_{\mathrm{interaction}}
#' \sim \mathrm{InvGamma}(a_\gamma, b_\gamma)} (a partial-pooling shrinkage
#' prior, shrinking towards the purely-additive model when the data don't
#' support an interaction). Defaults to \code{FALSE} (the original,
#' purely-additive mean).
#' @param gamma_proposal_window Proposal window (Gaussian random-walk SD)
#' for the interaction term; ignored if \code{include_interaction} is
#' \code{FALSE}.
#' @param a_gamma,b_gamma Shape/rate of the InvGamma hyperprior on the
#' interaction shrinkage variance; ignored if \code{include_interaction} is
#' \code{FALSE}.
#' @param batch_weight_prior One of \code{"global"} (the default),
#' \code{"partial_pooling"} or \code{"gp"}, controlling whether/how mixture
#' weights vary by batch:
#' \itemize{
#'   \item \code{"global"}: a single mixture weight vector shared by every
#'   batch - the original behaviour.
#'   \item \code{"partial_pooling"}: each batch gets its own weight vector,
#'   with the K-1 additive-log-ratio (ALR) coordinates of each batch's
#'   weights drawn exchangeably around a shared, estimated population mean
#'   and variance - no assumed order, distance or covariance structure
#'   between batches at all (the covariance is simply unknown/unstructured).
#'   Use this when batches' proportions are expected to vary but there is no
#'   known notion of which batches should be more similar to which.
#'   \item \code{"gp"}: each batch gets its own weight vector as above, but
#'   the ALR coordinates are instead linked across batches by a Gaussian
#'   process prior over \code{batch_coordinates} - for batches collected
#'   over a genuine, known ordering in time or space, where nearby batches
#'   are expected to be more similar than distant ones. Requires more
#'   assumptions than \code{"partial_pooling"} (an ordering/distance and a
#'   length scale), so prefer \code{"partial_pooling"} unless that ordering
#'   is actually known and relevant.
#' }
#' @param batch_coordinates A numeric vector of 1-D coordinates for the
#' batches (e.g. collection order or time), recycled/matched to the sorted
#' unique values of \code{batch_vec}. Defaults to \code{NULL}, i.e. batch
#' order 0, 1, ..., B - 1. Only used if \code{batch_weight_prior} is
#' \code{"gp"}.
#' @param gp_tau2,gp_length_scale Marginal variance and length scale of the
#' squared-exponential Gaussian process kernel over \code{batch_coordinates};
#' only used if \code{batch_weight_prior} is \code{"gp"}.
#' @param eta_proposal_window Proposal window for the batch-weight
#' Metropolis-Hastings update; used if \code{batch_weight_prior} is
#' \code{"partial_pooling"} or \code{"gp"}.
#' @param sample_gp_hyperparameters Logical; if \code{TRUE}, \code{gp_tau2}
#' and \code{gp_length_scale} are themselves updated by Metropolis-Hastings
#' rather than held fixed at the values passed in; only used if
#' \code{batch_weight_prior} is \code{"gp"}.
#' @param gp_hyperparameter_proposal_window Proposal window for the GP
#' hyperparameter update; only used if \code{batch_weight_prior} is
#' \code{"gp"} and \code{sample_gp_hyperparameters} is \code{TRUE}.
#' @param pp_tau2_shape,pp_tau2_rate Shape/rate of the InvGamma hyperprior
#' on each ALR coordinate's population variance (how much pooling there is
#' across batches: small values of the resulting tau2 pull batches strongly
#' towards their shared mean, large values let them vary close to
#' independently); only used if \code{batch_weight_prior} is
#' \code{"partial_pooling"}.
#' @param pp_mu_prior_sd Prior standard deviation for each ALR coordinate's
#' population mean (the shared value batches are pooled towards); only used
#' if \code{batch_weight_prior} is \code{"partial_pooling"}.
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
#' R <- 1000
#' thin <- 50
#'
#' # MCMC samples and BIC vector
#' samples <- batchSemiSupervisedMixtureModel(
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
#' initial_class_means <- matrix(c(1, 1, 3, 4), nrow = 2)
#' initial_class_covariance <- array(c(1, 0, 0, 1, 1, 0, 0, 1),
#'   dim = c(2, 2, 2)
#' )
#'
#' # We can use values from a previous chain
#' initial_batch_shift <- samples$batch_shift[, , R / thin]
#' initial_batch_scale <- matrix(
#'   c(1.2, 1.3, 1.7, 1.1, 1.4, 1.3, 1.2, 1.2, 1.1, 2.0),
#'   nrow = 2
#' )
#'
#' samples <- batchSemiSupervisedMixtureModel(X,
#'   R,
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
                                            R,
                                            thin,
                                            initial_labels,
                                            fixed,
                                            batch_vec,
                                            type,
                                            K_max = length(unique(initial_labels)),
                                            alpha = NULL,
                                            concentration = NULL,
                                            mu_proposal_window = 0.5**2,
                                            cov_proposal_window = 0.002,
                                            r_proposal_window = 0.1,
                                            sigma_proposal_window = 0.01,
                                            m_proposal_window = 0.3**2,
                                            S_proposal_window = 0.01,
                                            t_df_proposal_window = 0.015,
                                            m_scale = NULL,
                                            rho = 3.0,
                                            theta = 1.0,
                                            initial_class_means = NULL,
                                            initial_class_covariance = NULL,
                                            initial_batch_shift = NULL,
                                            initial_batch_scale = NULL,
                                            initial_class_df = NULL,
                                            verbose = TRUE,
                                            eta = 1.0,
                                            column_type = NULL,
                                            censor_code = NULL,
                                            auto_tune = TRUE,
                                            n_burn = NULL,
                                            include_interaction = FALSE,
                                            gamma_proposal_window = 0.1,
                                            a_gamma = 2.0,
                                            b_gamma = 1.0,
                                            batch_weight_prior = c("global", "partial_pooling", "gp"),
                                            batch_coordinates = NULL,
                                            gp_tau2 = 1.0,
                                            gp_length_scale = 1.0,
                                            eta_proposal_window = 0.1,
                                            sample_gp_hyperparameters = FALSE,
                                            gp_hyperparameter_proposal_window = 0.1,
                                            pp_tau2_shape = 2.0,
                                            pp_tau2_rate = 1.0,
                                            pp_mu_prior_sd = 10.0) {
  if (!is.matrix(X)) {
    stop("X is not a matrix. Data should be in matrix format.")
  }

  if (length(batch_vec) != nrow(X)) {
    stop("The number of rows in X and the number of batch labels are not equal.")
  }

  if (R < thin) {
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
    n_burn <- floor(R / 2)
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
      R,
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
      R,
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
      R,
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
      R,
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
  mcmc_output$R <- R
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
