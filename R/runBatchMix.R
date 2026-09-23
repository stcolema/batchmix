#' @title Run a batch mixture model (single chain)
#' @description \strong{The single canonical reference for every argument
#' shared across this package's model-fitting functions} -
#' \code{\link{fitBatchMix}} (the main, multi-chain entry point most users
#' should call) and \code{\link{batchSemiSupervisedMixtureModel}} (the
#' low-level engine underlying both, rarely called directly) document their
#' shared arguments by inheriting the descriptions below (\code{@@inheritParams
#' runBatchMix}) rather than repeating them, so there is exactly one place to
#' read about e.g. \code{control} - here. Runs a single MCMC chain
#' for a Bayesian mixture model which models both batch effects and
#' class/cluster structure (unsupervised if \code{fixed} is all 0 or not
#' given, semi-supervised otherwise). Every Metropolis-Hastings proposal
#' window (bundled into \code{control}, see \code{\link{batchmixControl}})
#' is auto-tuned by default (\code{control$auto_tune = TRUE}, via
#' Robbins-Monro diminishing adaptation over the first \code{control$n_burn}
#' iterations, frozen thereafter) - manually tuning proposal windows is not
#' required for typical use, and \code{control}'s fields are only starting
#' values for that adaptation, not fixed settings you need to get right
#' yourself.
#' @param X Data to cluster as a matrix with the items to cluster held in
#' rows. Missing entries (\code{NA}/\code{NaN}) are supported for every
#' \code{type}: each is modelled via proper per-sweep Gibbs data
#' augmentation (the missing entries of a row are redrawn from their full
#' conditional, given the other entries of that row and the current
#' cluster/batch parameters, every MCMC iteration - never imputed once and
#' held fixed, and never written back to \code{X} itself). This assumes the
#' missingness is at most Missing At Random given the modelled cluster/batch
#' structure and the other observed columns (MCAR is a special case) - not
#' Missing Not At Random (e.g. values missing because they are large/small),
#' which none of these sampler types attempt to model. The posterior draws
#' of the completed data are returned separately as
#' \code{mcmc_output$latent_data} (an N x P x saved-iterations array); the
#' original \code{X} you pass in is never modified. \code{type = "MVN_MIXED"}
#' additionally supports binary/probit columns and (non-missing) left/right
#' censored entries via \code{column_type}/\code{censor_code} - see the
#' 'probit_missing_censored' vignette.
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
#' \code{vignette("covariance_models", package = "batchmix")}), or
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
#' @param alpha The (symmetric) concentration parameter for the component
#' weights \code{w} in the model, used as \code{rep(alpha, K_max)} to build
#' \code{concentration}; only used if \code{concentration} (only available on
#' \code{\link{batchSemiSupervisedMixtureModel}} directly) is not given. Also
#' passed to \code{\link{generateInitialLabels}}'s stick-breaking prior
#' \emph{only} to draw a starting clustering when \code{initial_labels} is
#' not supplied - a different (GEM/Dirichlet-process) prior to the one above,
#' but since it only chooses where the chain starts and not the model being
#' fit, this does not bias the posterior; see
#' \code{\link{generateInitialLabels}}.
#' @param control A \code{\link{batchmixControl}} object bundling every
#' Metropolis-Hastings proposal window and the auto-tuning schedule
#' (\code{auto_tune}/\code{n_burn}) - the sampler-tuning knobs that rarely
#' need attention, as opposed to the prior/model-structure hyperparameters
#' below (\code{rho}, \code{theta}, \code{eta}, ...), which stay as
#' ordinary named arguments since they change what is being fitted, not how
#' hard the sampler works to fit it. See \code{?batchmixControl} for every
#' field and its default.
#' @param auto_tune,n_burn,mu_proposal_window,cov_proposal_window,r_proposal_window,sigma_proposal_window,m_proposal_window,S_proposal_window,t_df_proposal_window,gamma_proposal_window,eta_proposal_window,gp_hyperparameter_proposal_window
#' \strong{Deprecated}: pass these inside \code{control =
#' batchmixControl(...)} instead (e.g. \code{control =
#' batchmixControl(auto_tune = FALSE)} rather than \code{auto_tune =
#' FALSE}). Still work this release (with a warning); see
#' \code{?batchmixControl} for what each one does. If both \code{control}
#' and one of these are supplied, \code{control} wins and the deprecated
#' argument is ignored (with a warning).
#' @param verbose Logical indicating if warnings about proposal windows should
#' be printed.
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
                        control = batchmixControl(),
                        auto_tune = TRUE,
                        n_burn = NULL,
                        # -- proposal windows (deprecated - use `control` instead; only
                        # matter if auto_tune = FALSE) --
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

  deprecated_control_args <- list()
  if (!missing(auto_tune)) deprecated_control_args$auto_tune <- auto_tune
  if (!missing(n_burn)) deprecated_control_args$n_burn <- n_burn
  if (!missing(mu_proposal_window)) deprecated_control_args$mu_proposal_window <- mu_proposal_window
  if (!missing(cov_proposal_window)) deprecated_control_args$cov_proposal_window <- cov_proposal_window
  if (!missing(r_proposal_window)) deprecated_control_args$r_proposal_window <- r_proposal_window
  if (!missing(sigma_proposal_window)) deprecated_control_args$sigma_proposal_window <- sigma_proposal_window
  if (!missing(m_proposal_window)) deprecated_control_args$m_proposal_window <- m_proposal_window
  if (!missing(S_proposal_window)) deprecated_control_args$S_proposal_window <- S_proposal_window
  if (!missing(t_df_proposal_window)) deprecated_control_args$t_df_proposal_window <- t_df_proposal_window
  if (!missing(gamma_proposal_window)) deprecated_control_args$gamma_proposal_window <- gamma_proposal_window
  if (!missing(eta_proposal_window)) deprecated_control_args$eta_proposal_window <- eta_proposal_window
  if (!missing(gp_hyperparameter_proposal_window)) {
    deprecated_control_args$gp_hyperparameter_proposal_window <- gp_hyperparameter_proposal_window
  }
  control <- .resolveControlArgs(missing(control), control, deprecated_control_args, "runBatchMix")

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
    control = control,
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
    include_interaction = include_interaction,
    a_gamma = a_gamma,
    b_gamma = b_gamma,
    batch_weight_prior = batch_weight_prior,
    batch_coordinates = batch_coordinates,
    gp_tau2 = gp_tau2,
    gp_length_scale = gp_length_scale,
    sample_gp_hyperparameters = sample_gp_hyperparameters,
    pp_tau2_shape = pp_tau2_shape,
    pp_tau2_rate = pp_tau2_rate,
    pp_mu_prior_sd = pp_mu_prior_sd
  )

  mcmc_out
}
