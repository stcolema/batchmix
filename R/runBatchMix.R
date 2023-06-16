#' @title Run Batch Mixture Model
#' @description Runs a MCMC chain for a Bayesian mixture model which models
#' both batch effects and class/cluster structure.
#' @param X Data to cluster as a matrix with the items to cluster held in rows.
#' @param initial_labels Initial clustering.
#' @param batch_vec Labels identifying which batch each item being clustered is
#' from.
#' @param R The number of iterations in the sampler.
#' @param thin The factor by which the samples generated are thinned, e.g. if
#' ``thin=50`` only every 50th sample is kept.
#' @param type Character indicating density type to use. One of 'MVN'
#' (multivariate normal distribution) or 'MVT' (multivariate t distribution).
#' @param K_max The number of components to include (the upper bound on the
#' number of clusters in each sample). Defaults to the number of unique labels
#' in ``initial_labels``.
#' @param fixed Which items are fixed in their initial label. If not given,
#' defaults to a vector of 0 meaning the model is run unsupervised.
#' @param alpha The concentration parameter for the stick-breaking prior and the
#' weights in the model.
#' @param mu_proposal_window The proposal window for the cluster mean proposal
#' kernel. Making this smaller will normally increase the acceptance rate for
#' the proposed values in the Metropolis-Hastings sampler. The proposal density
#' is a Gaussian distribution, the window is the variance.
#' @param cov_proposal_window The proposal window for the cluster covariance
#' proposal kernel. The proposal density is a Wishart distribution, this
#' argument is the reciprocal of the degree of freedom. It is recommended to
#' set this aiming for accpetance rates of greater than 0.5 for the covariance
#' matrices (e.g., between 2e-03 and 1e-04 is a good range to consider
#' initially). As the entire covariance matrix is sampled at once exploration is
#' difficult.
#' @param m_proposal_window The proposal window for the batch mean proposal
#'  kernel. The proposal density is a Gaussian distribution, the window is the
#' variance.
#' @param S_proposal_window The proposal window for the batch standard deviation
#'  proposal kernel. The proposal density is a Gamma distribution, this
#' argument is the reciprocal of the rate. Recommended range to initially
#' consider is 0.015 to 2e-03, though smaller values might be necessary
#' particularly in higher dimensional data.
#' @param t_df_proposal_window The proposal window for the degrees of freedom
#' for the multivariate t distribution (not used if type is not 'MVT'). The
#' proposal density is a Gamma distribution, this argument is the reciprocal of
#' the rate. If the data is more Gaussian than the degrees of freedom might have
#' high acceptance rates regardless of the value chosen.
#' @param m_scale The scale hyperparameter for the batch shift prior
#' distribution. This defines the scale of the batch effect upon the mean and
#' should be in (0, 1]. If `NULL`, this quantity is sampled rather then fixed.
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
#' type <- "MVN"
#'
#' # Sampling parameters
#' R <- 1000
#' thin <- 50
#'
#' # MCMC samples
#' mcmc_out <- runBatchMix(
#'   X,
#'   R,
#'   thin,
#'   batch_vec,
#'   type,
#'   initial_labels = labels,
#'   fixed = fixed
#' )
#'
#' # Given an initial value for the parameters
#' initial_class_means <- matrix(c(1, 1, 3, 4), nrow = 2)
#' initial_class_covariance <- array(c(1, 0, 0, 1, 1, 0, 0, 1),
#'   dim = c(2, 2, 2)
#' )
#'
#' # We can use values from a previous chain
#' initial_batch_shift <- mcmc_out$batch_shift[, , R / thin]
#' initial_batch_scale <- matrix(
#'   c(1.2, 1.3, 1.7, 1.1, 1.4, 1.3, 1.2, 1.2, 1.1, 2.0),
#'   nrow = 2
#' )
#'
#' mcmc_out <- runBatchMix(X,
#'   R,
#'   thin,
#'   batch_vec,
#'   type,
#'   initial_labels = labels,
#'   fixed = fixed,
#'   initial_class_means = initial_class_means,
#'   initial_class_covariance = initial_class_covariance,
#'   initial_batch_shift = initial_batch_shift,
#'   initial_batch_scale = initial_batch_scale
#' )
#'
runBatchMix <- function(X,
                        R,
                        thin,
                        batch_vec,
                        type,
                        K_max = NULL,
                        initial_labels = NULL,
                        fixed = NULL,
                        alpha = 1,
                        mu_proposal_window = 0.5**2,
                        cov_proposal_window = 0.002,
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
                        verbose = TRUE) {
  unsupervised <- is.null(fixed)
  no_initial_partition_given <- is.null(initial_labels)
  if (!is.matrix(X)) {
    stop("X is not a matrix. Data should be in matrix format.")
  }

  semisupervised_with_no_initial_labels <- (!unsupervised) && no_initial_partition_given

  if (semisupervised_with_no_initial_labels) {
    .err <- paste0(
      "If running a semi-supervised model (fixed has any entries ",
      "equal to 1), then an initial labelling must be given."
    )
    stop(.err)
  }

  N <- nrow(X)
  if (unsupervised) {
    fixed <- rep(0, N)
  }

  no_k_passed <- is.null(K_max)
  if (no_k_passed) {
    if (no_initial_partition_given) {
      K_max <- min(floor(N / 2), 50)
    } else {
      K_max <- length(unique(initial_labels))
    }
  }

  if (no_initial_partition_given) {
    initial_labels <- generateInitialLabels(alpha, K_max, fixed)
  }

  mcmc_out <- batchSemiSupervisedMixtureModel(X,
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

  mcmc_out
}
