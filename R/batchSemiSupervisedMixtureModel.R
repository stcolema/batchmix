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
#' (multivariate normal distribution) or 'MVT' (multivariate t distribution).
#' @param batch_specific_weights Allow each batch to have unique class weights.
#' If FALSE class weights are common across batches, if TRUE the class weights
#' are nested within batch and share a common hyperparameter. Defaults to TRUE.
#' @param K_max The number of components to include (the upper bound on the
#' number of clusters in each sample). Defaults to the number of unique labels
#' in ``initial_labels``.
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
#' samples <- batchSemiSupervisedMixtureModel(X,
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
                                            batch_specific_weights = TRUE,
                                            K_max = length(unique(initial_labels)),
                                            alpha = NULL,
                                            concentration = NULL,
                                            mu_proposal_window = 0.5**2,
                                            cov_proposal_window = 0.002,
                                            m_proposal_window = 0.3**2,
                                            S_proposal_window = 0.01,
                                            t_df_proposal_window = 0.015,
                                            m_scale = 0.1,
                                            rho = 11.0,
                                            theta = 5.0,
                                            initial_class_means = NULL,
                                            initial_class_covariance = NULL,
                                            initial_batch_shift = NULL,
                                            initial_batch_scale = NULL,
                                            initial_class_df = NULL,
                                            verbose = TRUE
) {
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
  if((! alpha_not_passed) && (! concentration_not_passed)) {
    stop("Only one of ``concentration`` or ``alpha`` should be passed.")
  }
  if(concentration_not_passed) {
    if (alpha_not_passed) {
      alpha = 1.0
    }
    concentration <- rep(alpha, K_max)
  }
  #   concentration <- table(initial_labels[fixed == 1]) / sum(fixed)
  # } else {
  #   concentration <- rep(1, K_max)
  # }
  
  # Check the proposal windows are all strictly positive
  checkProposalWindows(mu_proposal_window,
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
  
  P <- ncol(X)
  
  # Check if an initial value is passed for any of the parameters. Prepare the
  # parameters to be passed to C++.
  class_mean_passed <- !is.null(initial_class_means)
  class_covariance_passed <- !is.null(initial_class_covariance)
  batch_shift_passed <- !is.null(initial_batch_shift)
  batch_scale_passed <- !is.null(initial_batch_scale)
  class_df_passed <- !is.null(initial_class_df)
  
  initial_parameters <- prepareInitialParameters(initial_class_means,
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
  
  # Pull samples from the mixture model
  if(batch_specific_weights) {
    if (type == "MVN") {
      mcmc_output <- sampleSemisupervisedMVNVaryingWeights(
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
        batch_shift_passed,
        batch_scale_passed
      )
    }
    
    if (type == "MVT") {
      mcmc_output <- sampleSemisupervisedMVTVaryingWeights(
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
        batch_shift_passed,
        batch_scale_passed
      )
    }
    
    mcmc_output$weights <- t(apply(mcmc_output$weights, 3, function(x) c(x)))
    
  } else {
  if (type == "MVN") {
    mcmc_output <- sampleSemisupervisedMVN(
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
      batch_shift_passed,
      batch_scale_passed
    )
  }
  
  if (type == "MVT") {
    mcmc_output <- sampleSemisupervisedMVT(
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
      batch_shift_passed,
      batch_scale_passed
    )
  }
    
    mcmc_output$concentration <- matrix(concentration, nrow = 1)
  }
  
  if (!type %in% c("MVN", "MVT")) {
    stop("Type not recognised. Please use one of 'MVN' or 'MVT'.")
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
  mcmc_output$m_proposal_window <- m_proposal_window
  mcmc_output$S_proposal_window <- S_proposal_window
  mcmc_output$t_df_proposal_window <- t_df_proposal_window
  
  # Indicate if the model was semi-supervised or unsupervised
  mcmc_output$Semisupervised <- TRUE
  
  # Correct this if we were effectively unsupervised
  actually_unsupervised <- all(fixed == 0)
  if (actually_unsupervised) {
    mcmc_output$Semisupervised <- FALSE
  }
  
  # Record if component weights vary across batches
  mcmc_output$batch_specific_weights <- batch_specific_weights
  
  mcmc_output
}
