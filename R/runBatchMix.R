#' @title Run a batch mixture model (single chain)
#' @description \strong{The single canonical reference for every argument
#' shared across this package's model-fitting functions} -
#' \code{\link{fitBatchMix}} (the main, multi-chain entry point most users
#' should call) and \code{\link{batchSemiSupervisedMixtureModel}} (the
#' low-level engine underlying both, rarely called directly) document their
#' shared arguments by inheriting the descriptions below (\code{@@inheritParams
#' runBatchMix}) rather than repeating them, so there is exactly one place to
#' read about e.g. \code{mu_proposal_window} - here. Runs a single MCMC chain
#' for a Bayesian mixture model which models both batch effects and
#' class/cluster structure (unsupervised if \code{fixed} is all 0 or not
#' given, semi-supervised otherwise). Every Metropolis-Hastings proposal
#' window is auto-tuned by default (\code{auto_tune = TRUE}, via
#' Robbins-Monro diminishing adaptation over the first \code{n_burn}
#' iterations, frozen thereafter) - manually tuning \code{mu_proposal_window}
#' and friends is not required for typical use, and the proposal-window
#' arguments below are only starting values for that adaptation, not fixed
#' settings you need to get right yourself.
#' @param X Data to cluster as a matrix with the items to cluster held in
#' rows. Missing entries (\code{NA}/\code{NaN}) are only supported when
#' \code{type = "MVN_MIXED"}, which models them via proper data augmentation
#' (see \code{column_type}/\code{censor_code} below and the
#' 'probit_missing_censored' vignette); for every other \code{type}, X must
#' be complete - the empirical-Bayes prior setup for those samplers uses
#' \code{mean(X)}/\code{cov(X)} directly and does not tolerate missing
#' values.
#' @param n_iter The number of iterations in the sampler.
#' @param thin The factor by which the samples generated are thinned, e.g. if
#' ``thin=50`` only every 50th sample is kept.
#' @param batch_vec Labels identifying which batch each item being clustered is
#' from.
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
#' in ``initial_labels`` if given, else \code{min(floor(N / 2), 50)}.
#' @param initial_labels Initial clustering. If not given, defaults to a
#' k-means partition of \code{X} for a fully unsupervised fit, or a random
#' draw from the stick-breaking prior otherwise - see
#' \code{\link{generateInitialLabels}}.
#' @param fixed Which items are fixed in their initial label. If not given,
#' defaults to a vector of 0 meaning the model is run unsupervised.
#' @param alpha The concentration parameter for the stick-breaking prior and
#' the weights in the model. Only used if \code{concentration} (only
#' available on \code{\link{batchSemiSupervisedMixtureModel}} directly) is
#' not given.
#' @param auto_tune Logical; if \code{TRUE} (the default), every proposal
#' window is adapted during the first \code{n_burn} iterations via
#' Robbins-Monro diminishing adaptation, instead of staying fixed at the
#' value passed in. Adaptation is frozen after \code{n_burn} iterations so
#' the post-burn-in chain retains the correct stationary distribution.
#' @param n_burn Number of iterations treated as burn-in for proposal-window
#' adaptation; ignored if \code{auto_tune} is \code{FALSE}. Defaults to half
#' of \code{n_iter}.
#' @param verbose Logical indicating if warnings about proposal windows should
#' be printed.
#' @param mu_proposal_window The proposal window for the cluster mean proposal
#' kernel. The proposal density is a Gaussian distribution, the window is the
#' variance. Making this smaller will normally increase the acceptance rate.
#' @param cov_proposal_window The proposal window for the cluster covariance
#' proposal kernel when \code{type} is 'MVN' or 'MVT'. The proposal density
#' is a Wishart distribution, this argument is the reciprocal of the degree
#' of freedom. It is recommended to aim for acceptance rates greater than
#' 0.5 (e.g. between 2e-03 and 1e-04 is a good range to consider initially) -
#' as the entire covariance matrix is sampled at once, exploration is
#' difficult.
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
#' argument is the reciprocal of the rate. Recommended range to initially
#' consider is 0.015 to 2e-03, though smaller values might be necessary
#' particularly in higher dimensional data.
#' @param t_df_proposal_window The proposal window for the degrees of freedom
#' for the multivariate t distribution (not used if type is not 'MVT'). The
#' proposal density is a Gamma distribution, this argument is the reciprocal of
#' the rate. If the data is more Gaussian than the degrees of freedom might have
#' high acceptance rates regardless of the value chosen.
#' @param gamma_proposal_window Proposal window (Gaussian random-walk SD)
#' for the interaction term; ignored if \code{include_interaction} is
#' \code{FALSE}.
#' @param eta_proposal_window Proposal window for the batch-weight
#' Metropolis-Hastings update; used if \code{batch_weight_prior} is
#' \code{"partial_pooling"} or \code{"gp"}.
#' @param gp_hyperparameter_proposal_window Proposal window for the GP
#' hyperparameter update; only used if \code{batch_weight_prior} is
#' \code{"gp"} and \code{sample_gp_hyperparameters} is \code{TRUE}.
#' @param m_scale The scale hyperparameter for the batch shift prior
#' distribution. This defines the scale of the batch effect upon the mean and
#' should be in (0, 1]. If `NULL`, this quantity is sampled rather then fixed.
#' @param rho The shape of the prior distribution for the batch scale.
#' @param theta The scale of the prior distribution for the batch scale.
#' @param eta The LKJ concentration parameter for the cluster correlation
#' prior; only used if \code{type} is 'MVN_LKJ' or 'MVN_MIXED'. eta = 1 is
#' uniform over the space of correlation matrices, eta > 1 shrinks
#' correlations towards 0.
#' @param a_gamma,b_gamma Shape/rate of the InvGamma hyperprior on the
#' interaction shrinkage variance; ignored if \code{include_interaction} is
#' \code{FALSE}.
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
#' @param include_interaction Logical; if \code{TRUE}, add a batch x cluster
#' interaction term to the mean, \eqn{\gamma_{k,b} \sim N(0,
#' \tau^2_{\mathrm{interaction}})}, with \eqn{\tau^2_{\mathrm{interaction}}
#' \sim \mathrm{InvGamma}(a_\gamma, b_\gamma)} (a partial-pooling shrinkage
#' prior, shrinking towards the purely-additive model when the data don't
#' support an interaction). Defaults to \code{FALSE} (the original,
#' purely-additive mean).
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
#' @param sample_gp_hyperparameters Logical; if \code{TRUE}, \code{gp_tau2}
#' and \code{gp_length_scale} are themselves updated by Metropolis-Hastings
#' rather than held fixed at the values passed in; only used if
#' \code{batch_weight_prior} is \code{"gp"}.
#' @param pp_tau2_shape,pp_tau2_rate Shape/rate of the InvGamma hyperprior
#' on each ALR coordinate's population variance (how much pooling there is
#' across batches: small values of the resulting tau2 pull batches strongly
#' towards their shared mean, large values let them vary close to
#' independently); only used if \code{batch_weight_prior} is
#' \code{"partial_pooling"}.
#' @param pp_mu_prior_sd Prior standard deviation for each ALR coordinate's
#' population mean (the shared value batches are pooled towards); only used
#' if \code{batch_weight_prior} is \code{"partial_pooling"}.
#' @param column_type Only used if \code{type} is 'MVN_MIXED': a P-vector, 0
#' for a continuous column, 1 for a binary column observed via a probit
#' link. Defaults to all-continuous.
#' @param censor_code Only used if \code{type} is 'MVN_MIXED': an N x P
#' matrix, meaningful only for continuous columns: 0 = not censored,
#' 1 = left-censored (true value below the recorded X entry), 2 =
#' right-censored (true value above the recorded X entry). Defaults to no
#' censoring. Missing entries are indicated by \code{NA}/\code{NaN} in
#' \code{X} directly, for both continuous and binary columns.
#' @param ... Accepts only the deprecated \code{R} argument (renamed to
#' \code{n_iter}; still works, with a warning, for this release, when
#' passed positionally in its original slot or when every other argument is
#' also named - naming \code{R} while leaving later arguments positional is
#' not supported, since there is no way to bind two names to the same
#' argument slot). Anything else is an "unused argument" error.
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
#' n_iter <- 1000
#' thin <- 50
#'
#' # MCMC samples
#' mcmc_out <- runBatchMix(
#'   X,
#'   n_iter,
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
#' initial_batch_shift <- mcmc_out$batch_shift[, , n_iter / thin]
#' initial_batch_scale <- matrix(
#'   c(1.2, 1.3, 1.7, 1.1, 1.4, 1.3, 1.2, 1.2, 1.1, 2.0),
#'   nrow = 2
#' )
#'
#' mcmc_out <- runBatchMix(X,
#'   n_iter,
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
                        n_iter,
                        thin,
                        batch_vec,
                        type,
                        # -- problem specification --
                        K_max = NULL,
                        initial_labels = NULL,
                        fixed = NULL,
                        alpha = 1,
                        # -- MCMC control (auto-tuning is on by default - see Description) --
                        auto_tune = TRUE,
                        n_burn = NULL,
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
                        verbose = TRUE,
                        ...) {
  n_iter <- .resolveDeprecatedNIter(
    missing(n_iter), if (missing(n_iter)) NULL else n_iter, list(...), "runBatchMix"
  )

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
    initial_labels <- generateInitialLabels(alpha, K_max, fixed, X = X)
  }

  mcmc_out <- batchSemiSupervisedMixtureModel(X,
    n_iter,
    thin,
    initial_labels,
    fixed,
    batch_vec,
    type,
    K_max = K_max,
    alpha = alpha,
    mu_proposal_window = mu_proposal_window,
    cov_proposal_window = cov_proposal_window,
    r_proposal_window = r_proposal_window,
    sigma_proposal_window = sigma_proposal_window,
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
    verbose = verbose,
    eta = eta,
    column_type = column_type,
    censor_code = censor_code,
    auto_tune = auto_tune,
    n_burn = n_burn,
    include_interaction = include_interaction,
    gamma_proposal_window = gamma_proposal_window,
    a_gamma = a_gamma,
    b_gamma = b_gamma,
    batch_weight_prior = batch_weight_prior,
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

  mcmc_out
}
