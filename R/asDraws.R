# Conversion of fits to the draws formats of the posterior package, so that the
# summaries, diagnostics and plots of the Stan workflow (posterior, bayesplot)
# apply to batchmix output. The methods are registered when posterior is loaded.

# Internal: the indices of the stored samples of one chain that are retained.
# This is the rule of assessConvergence() and calcBICM(): a chain already
# trimmed by processMCMCChain() (it carries `pred`, or has fewer draws than the
# sampler saved) keeps everything; otherwise the burn in is
# `burn` iterations (default half of `n_iter`), whatever the thinning.
.batchmixRetainedDraws <- function(chain, burn = NULL) {
  n_saved <- nrow(chain$samples)
  n_iter <- chain$n_iter
  thin <- chain$thin
  # processMCMCChain() adds `pred` and removes the burn in
  if (!is.null(chain$pred) || n_saved < floor(n_iter / thin)) {
    return(seq_len(n_saved))
  }
  if (is.null(burn)) burn <- floor(n_iter / 2)
  first <- max(1L, ceiling((burn + thin) / thin))
  if (first > n_saved) {
    stop("A burn in of ", burn, " iterations leaves no saved samples (", n_saved,
      " were saved, thin = ", thin, ").",
      call. = FALSE
    )
  }
  seq.int(first, n_saved)
}

# Internal: the scalar quantities that do not depend on how the clusters are
# labelled (or, optionally, label-dependent ones the user asks for), as a draws
# x variables matrix for the stored samples `index` of one chain.
.batchmixDrawsMatrix <- function(x, index, batch = TRUE, weights = FALSE, allocations = FALSE) {
  n <- length(index)
  labels <- x$samples[index, , drop = FALSE]
  cols <- list()
  cols[["log_lik_complete"]] <- as.numeric(x$complete_likelihood)[index]
  if (!is.null(x$observed_likelihood)) {
    cols[["log_lik_observed"]] <- as.numeric(x$observed_likelihood)[index]
  }
  cols[["bic"]] <- as.numeric(x$BIC)[index]
  cols[["occupied"]] <- apply(labels, 1, function(z) length(unique(z)))

  # Hyperparameters, only when the sampler updated them (otherwise they are constants)
  if (isTRUE(x$sample_m_scale)) cols[["lambda_2"]] <- as.numeric(x$lambda_2)[index]
  if (isTRUE(x$sample_s_scale)) cols[["rho"]] <- as.numeric(x$rho)[index]
  if (identical(x$batch_weight_prior, "gp") && isTRUE(x$sample_gp_hyperparameters)) {
    cols[["gp_tau2"]] <- as.numeric(x$gp_tau2)[index]
    # the non-stationary kernels have no length scale
    if (identical(x$gp_kernel, "matern32")) {
      cols[["gp_length_scale"]] <- as.numeric(x$gp_length_scale)[index]
    }
  }

  if (batch) {
    # Batches are observed, so these mean the same thing in every draw
    P <- dim(x$batch_shift)[1]
    B <- dim(x$batch_shift)[2]
    for (nm in c("batch_shift", "batch_scale")) {
      arr <- x[[nm]][, , index, drop = FALSE]
      for (b in seq_len(B)) {
        for (p in seq_len(P)) {
          cols[[sprintf("%s[%d,%d]", nm, p, b)]] <- arr[p, b, ]
        }
      }
    }
  }

  if (weights) {
    # Normalised weights sorted within each draw, so that they do not depend on the labelling
    sort_rows <- function(w) {
      w <- w / rowSums(w)
      matrix(t(apply(w, 1, sort, decreasing = TRUE)), nrow = nrow(w))
    }
    if (is.null(x$batch_weight_prior) || x$batch_weight_prior == "global") {
      w <- sort_rows(matrix(x$weights[index, ], nrow = n))
      for (k in seq_len(ncol(w))) cols[[sprintf("sorted_weight[%d]", k)]] <- w[, k]
    } else {
      # one set of weights per batch: draws x K from the B x K x draws array
      for (b in seq_len(dim(x$w_batch)[1])) {
        w <- sort_rows(matrix(x$w_batch[b, , index], nrow = n, byrow = TRUE))
        for (k in seq_len(ncol(w))) cols[[sprintf("sorted_weight[%d,%d]", b, k)]] <- w[, k]
      }
    }
  }

  if (allocations) {
    # Labelled-dependent unless the items are observed; numbered from 1
    z <- labels + 1
    for (i in seq_len(ncol(z))) cols[[sprintf("allocation[%d]", i)]] <- z[, i]
  }
  do.call(cbind, cols)
}

.batchmixDrawsArray <- function(chains, burn = NULL, batch = TRUE, weights = FALSE, allocations = FALSE) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("This needs the 'posterior' package.", call. = FALSE)
  }
  mats <- lapply(chains, function(ch) {
    .batchmixDrawsMatrix(ch, .batchmixRetainedDraws(ch, burn), batch, weights, allocations)
  })
  n_draws <- vapply(mats, nrow, integer(1))
  if (length(unique(n_draws)) != 1) {
    stop("The chains retain different numbers of draws.", call. = FALSE)
  }
  if (length(unique(lapply(mats, colnames))) != 1) {
    stop("The chains were fitted with different settings, so they have different variables.", call. = FALSE)
  }
  arr <- array(
    unlist(mats, use.names = FALSE),
    dim = c(n_draws[1], ncol(mats[[1]]), length(mats)),
    dimnames = list(iteration = NULL, variable = colnames(mats[[1]]), chain = NULL)
  )
  posterior::as_draws_array(aperm(arr, c(1, 3, 2)))
}

.batchmixChains <- function(x) {
  if (inherits(x, "batchmix_fit_list")) unclass(x)[seq_along(x)] else list(x)
}

#' @title Convert a batchmix fit to the draws formats of the posterior package
#' @description Methods of \code{posterior::as_draws()} and its relatives for
#' the output of \code{\link{runBatchMix}},
#' \code{\link{batchSemiSupervisedMixtureModel}}, \code{\link{fitBatchMix}} and
#' \code{\link{processMCMCChains}}. With the result, the Stan workflow carries
#' over: \code{posterior::summarise_draws()} (mean, quantiles, Rhat, bulk and
#' tail ESS), \code{posterior::subset_draws()}, \code{posterior::resample_draws()},
#' and the \pkg{bayesplot} functions that take draws.
#'
#' A mixture model has no parameters that mean the same thing in every draw
#' (clusters can swap labels), so the variables are restricted to quantities
#' that do not depend on the labelling:
#' \describe{
#'   \item{\code{log_lik_complete}, \code{log_lik_observed}, \code{bic}}{The
#'   log-likelihoods and BIC recorded by the sampler (see
#'   \code{\link{assessConvergence}}). The log-likelihoods are not \code{lp__}:
#'   they omit the priors.}
#'   \item{\code{occupied}}{Number of occupied clusters.}
#'   \item{\code{batch_shift[p,b]}, \code{batch_scale[p,b]}}{The shift and scale
#'   of batch \code{b} in feature \code{p}. Batches are observed, so these keep
#'   their meaning across draws. Dropped with \code{batch = FALSE}.}
#'   \item{\code{lambda_2}, \code{rho}}{Hyperparameters of the batch shift and
#'   batch scale priors, present only if the sampler updated them
#'   (\code{sample_m_scale}, \code{sample_s_scale}).}
#'   \item{\code{gp_tau2}, \code{gp_length_scale}}{Gaussian process
#'   hyperparameters, present only for \code{batch_weight_prior = "gp"} with
#'   \code{sample_gp_hyperparameters = TRUE}; the length scale only for the
#'   \code{"matern32"} kernel, since the other kernels have none.}
#'   \item{\code{sorted_weight[k]} or \code{sorted_weight[b,k]}}{With
#'   \code{weights = TRUE}: the normalised mixture weights sorted from largest
#'   to smallest within each draw, for each batch under batch-specific weights.}
#'   \item{\code{allocation[i]}}{With \code{allocations = TRUE}: the cluster
#'   (from 1) of item \code{i}. These depend on the labelling, so their Rhat and
#'   ESS are meaningful only for items whose labels are observed.}
#' }
#' Cluster means, covariances and weights are left out for the same reason; the
#' interaction terms and per-cluster hyperparameters (\code{pp_mu}, \code{pp_tau2})
#' are cluster-indexed too.
#'
#' The burn in follows \code{\link{assessConvergence}}: the first \code{burn}
#' iterations are dropped, so the diagnostics of \pkg{posterior} agree with it.
#' Chains that have been through \code{\link{processMCMCChain}} are already trimmed
#' and keep every draw.
#' @param x Output of \code{\link{runBatchMix}} or
#' \code{\link{batchSemiSupervisedMixtureModel}} (\code{batchmix_fit}), or of
#' \code{\link{fitBatchMix}} / \code{\link{processMCMCChains}}
#' (\code{batchmix_fit_list}).
#' @param burn Number of iterations treated as burn in. The default is half of
#' \code{n_iter}. Ignored for chains that have been through
#' \code{\link{processMCMCChain}}.
#' @param batch Include the batch shifts and scales.
#' @param weights Include the sorted normalised mixture weights.
#' @param allocations Include the item allocations (large: one variable per item).
#' @param ... Unused.
#' @return A draws object of the requested class: \code{as_draws()} and
#' \code{as_draws_array()} give a \code{draws_array} (iterations x chains x
#' variables), the others the corresponding format.
#' @seealso \code{\link{assessConvergence}}, \code{\link{rankNormalizedRhat}}
#' @name as_draws.batchmix_fit
#' @examples
#' if (requireNamespace("posterior", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- rbind(matrix(rnorm(80, 0, 1), ncol = 2), matrix(rnorm(80, 3, 1), ncol = 2))
#'   batch_vec <- rep(1:2, length.out = nrow(X))
#'   fit <- fitBatchMix(X,
#'     n_chains = 2, n_iter = 400, thin = 5, batch_vec = batch_vec,
#'     type = "MVN", K_max = 4, verbose = FALSE
#'   )
#'   draws <- posterior::as_draws(fit)
#'   posterior::summarise_draws(draws)
#' }
NULL

#' @rdname as_draws.batchmix_fit
#' @exportS3Method posterior::as_draws
as_draws.batchmix_fit <- function(x, burn = NULL, batch = TRUE, weights = FALSE, allocations = FALSE, ...) {
  .batchmixDrawsArray(.batchmixChains(x), burn, batch, weights, allocations)
}

#' @rdname as_draws.batchmix_fit
#' @exportS3Method posterior::as_draws
as_draws.batchmix_fit_list <- as_draws.batchmix_fit

#' @rdname as_draws.batchmix_fit
#' @exportS3Method posterior::as_draws_array
as_draws_array.batchmix_fit <- as_draws.batchmix_fit

#' @rdname as_draws.batchmix_fit
#' @exportS3Method posterior::as_draws_array
as_draws_array.batchmix_fit_list <- as_draws.batchmix_fit

#' @rdname as_draws.batchmix_fit
#' @exportS3Method posterior::as_draws_df
as_draws_df.batchmix_fit <- function(x, ...) posterior::as_draws_df(as_draws.batchmix_fit(x, ...))

#' @rdname as_draws.batchmix_fit
#' @exportS3Method posterior::as_draws_df
as_draws_df.batchmix_fit_list <- as_draws_df.batchmix_fit

#' @rdname as_draws.batchmix_fit
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.batchmix_fit <- function(x, ...) posterior::as_draws_matrix(as_draws.batchmix_fit(x, ...))

#' @rdname as_draws.batchmix_fit
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.batchmix_fit_list <- as_draws_matrix.batchmix_fit

#' @rdname as_draws.batchmix_fit
#' @exportS3Method posterior::as_draws_list
as_draws_list.batchmix_fit <- function(x, ...) posterior::as_draws_list(as_draws.batchmix_fit(x, ...))

#' @rdname as_draws.batchmix_fit
#' @exportS3Method posterior::as_draws_list
as_draws_list.batchmix_fit_list <- as_draws_list.batchmix_fit

#' @rdname as_draws.batchmix_fit
#' @exportS3Method posterior::as_draws_rvars
as_draws_rvars.batchmix_fit <- function(x, ...) posterior::as_draws_rvars(as_draws.batchmix_fit(x, ...))

#' @rdname as_draws.batchmix_fit
#' @exportS3Method posterior::as_draws_rvars
as_draws_rvars.batchmix_fit_list <- as_draws_rvars.batchmix_fit
