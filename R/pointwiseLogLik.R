#' @title Pointwise log-likelihood of the fitted items
#' @description The log-likelihood of each item at each retained draw, with the
#' item's cluster allocation summed out. The result is the draws x items matrix
#' that \code{loo::waic()} and \code{loo::loo()} take (Vehtari, Gelman and
#' Gabry, 2017); for several chains it carries the chain of each row in
#' \code{attr(., "chain_id")}, for \code{loo::relative_eff()}.
#'
#' For an item \eqn{n} in batch \eqn{b} the value at a draw is
#' \deqn{\log \sum_k w_{bk}\, f(x_n \mid \mu_{kb}, \Sigma_{kb}),}
#' where \eqn{w_b} is the weight vector of batch \eqn{b} (the global weights
#' when \code{batch_weight_prior = "global"}) and \eqn{f} is the multivariate
#' normal density (\code{type} \code{"MVN"}, \code{"MVN_LKJ"}) or multivariate
#' t density (\code{"MVT"}, with the cluster's degrees of freedom and
#' \eqn{\Sigma_{kb}} as scale matrix) of the cluster's mean and covariance after
#' the batch shift, the batch scale and, if used, the batch-by-cluster
#' interaction have been applied. Because the allocation is summed out, the
#' value depends on the parameters only, not on any item-specific latent
#' variable, which is the form that information criteria and leave-one-out
#' cross-validation need. It also does not depend on how the clusters are
#' labelled.
#'
#' The values are recomputed from the stored draws, so no refit is needed. They
#' need the data and the batch labels, which a fit does not store. Only the
#' multivariate normal and t types are supported: \code{type = "MVN_MIXED"} has
#' binary and censored columns whose likelihood involves probit integrals, which
#' are not computed here.
#'
#' An entry of \code{X} that is \code{NA} is integrated out: an item contributes
#' the density of its observed coordinates (the marginal of the normal or t
#' distribution is the same distribution on those coordinates), and an item with
#' no observed coordinate contributes 0.
#'
#' Items whose label is observed (semi-supervised fits, \code{fixed}) contribute
#' the joint density of their data and of that label,
#' \eqn{\log w_{b y_n} f(x_n \mid \cdot)}, as in the observed-data
#' log-likelihood the sampler records. They are not comparable with items whose
#' label is summed out; for leave-one-out comparison of models with different
#' sets of labelled items, use the unlabelled items only.
#'
#' \strong{What leave-one-out means here.} Leaving one item out keeps its batch
#' in the data, so the batch shifts, scales and weights are estimated from the
#' other items of the batch. This measures how well the model predicts a new
#' item from a batch it has seen. It does not measure prediction for an unseen
#' batch; see \code{\link{predictNewBatch}} for that. Pareto-\eqn{k} diagnostics
#' from \pkg{loo} should be checked: an item that is the only member of a
#' cluster, or of a batch, is influential and its importance sampling
#' approximation may be unreliable.
#' @param object Output of \code{\link{runBatchMix}} or
#' \code{\link{batchSemiSupervisedMixtureModel}} (\code{batchmix_fit}), or of
#' \code{\link{fitBatchMix}} (\code{batchmix_fit_list}), before or after
#' \code{\link{processMCMCChain}}.
#' @param X The data the model was fitted to: an N x P numeric matrix, possibly
#' with \code{NA} entries.
#' @param batch_vec The batch labels used in the fit, of length N.
#' @param burn Number of iterations treated as burn in. The default is half of
#' \code{n_iter}, as in \code{\link{assessConvergence}} and
#' \code{\link{as_draws.batchmix_fit}}. Ignored for chains that have been
#' through \code{\link{processMCMCChain}}.
#' @param fixed For a semi-supervised fit, the binary vector (length N) of the
#' items whose label was observed, as passed to the fit. Required if the fit
#' was semi-supervised, and not used otherwise.
#' @return A matrix with a row per retained draw (all chains stacked) and a
#' column per item, named by \code{rownames(X)} if it has any. The chain of each
#' row is in \code{attr(., "chain_id")}.
#' @references Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian
#' model evaluation using leave-one-out cross-validation and WAIC.
#' \emph{Statistics and Computing}, 27, 1413-1432.
#' @seealso \code{\link{predictNewBatch}}, \code{\link{as_draws.batchmix_fit}}
#' @export
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(80, 0, 1), ncol = 2), matrix(rnorm(80, 3, 1), ncol = 2))
#' batch_vec <- rep(1:2, length.out = nrow(X))
#' fit <- fitBatchMix(X,
#'   n_chains = 2, n_iter = 200, thin = 5, batch_vec = batch_vec,
#'   type = "MVN", K_max = 4, verbose = FALSE
#' )
#' ll <- pointwiseLogLik(fit, X, batch_vec)
#' dim(ll)
#'
#' if (requireNamespace("loo", quietly = TRUE)) {
#'   r_eff <- loo::relative_eff(exp(ll), chain_id = attr(ll, "chain_id"))
#'   loo::loo(ll, r_eff = r_eff)
#' }
pointwiseLogLik <- function(object, X, batch_vec, burn = NULL, fixed = NULL) {
  chains <- .batchmixChains(object)
  if (!is.list(chains[[1]]) || is.null(chains[[1]]$mean_sum)) {
    stop("`object` must be the output of runBatchMix(), batchSemiSupervisedMixtureModel() or fitBatchMix().",
      call. = FALSE
    )
  }
  first <- chains[[1]]
  if (identical(first$type, "MVN_MIXED")) {
    stop("pointwiseLogLik() does not support type = 'MVN_MIXED': the binary and censored columns ",
      "need probit integrals that are not computed here.",
      call. = FALSE
    )
  }
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X) != first$N || ncol(X) != first$P) {
    stop("`X` must be a numeric ", first$N, " x ", first$P, " matrix, as in the fit.", call. = FALSE)
  }
  if (length(batch_vec) != first$N) {
    stop("`batch_vec` must have one entry per row of `X`.", call. = FALSE)
  }
  # The recoding used when fitting and in simulatePosteriorPredictive()
  if (!any(batch_vec == 0)) {
    batch_vec <- as.numeric(as.factor(batch_vec)) - 1
  }
  if (length(unique(batch_vec)) != first$B) {
    stop("`batch_vec` has ", length(unique(batch_vec)), " batches, but the fit has ", first$B, ".", call. = FALSE)
  }
  if (isTRUE(first$Semisupervised)) {
    if (is.null(fixed)) {
      stop("This fit was semi-supervised: give `fixed`, the items whose label was observed.", call. = FALSE)
    }
    if (length(fixed) != first$N) stop("`fixed` must have one entry per row of `X`.", call. = FALSE)
    fixed <- as.logical(fixed)
  } else {
    fixed <- rep(FALSE, first$N)
  }

  groups <- .batchmixItemGroups(X, batch_vec)
  n_regularised <- 0L
  out <- lapply(chains, function(ch) {
    res <- .batchmixPointwiseChain(ch, X, groups, fixed, .batchmixRetainedDraws(ch, burn))
    n_regularised <<- n_regularised + res$n_regularised
    res$ll
  })
  if (n_regularised > 0) {
    warning(n_regularised, " covariance matrices were not numerically positive definite and were regularised.",
      call. = FALSE
    )
  }
  chain_id <- rep(seq_along(out), vapply(out, nrow, integer(1)))
  ll <- do.call(rbind, out)
  colnames(ll) <- rownames(X)
  attr(ll, "chain_id") <- chain_id
  ll
}

# Items grouped by batch and pattern of observed coordinates, so that each group
# shares one set of component densities
.batchmixItemGroups <- function(X, batch_vec) {
  observed <- !is.na(X)
  key <- paste(batch_vec, apply(observed, 1, function(r) paste(as.integer(r), collapse = "")))
  lapply(split(seq_len(nrow(X)), key), function(items) {
    list(items = items, batch = batch_vec[items[1]], observed = which(observed[items[1], ]))
  })
}

# Cholesky factor, regularised if the matrix is not numerically positive
# definite (a Metropolis proposal can leave a stored covariance on the boundary).
# `flag` counts the regularisations.
.batchmixCholSafe <- function(S, flag) {
  ch <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(ch)) {
    eig <- eigen((S + t(S)) / 2, symmetric = TRUE)
    vals <- pmax(eig$values, .Machine$double.eps * max(abs(eig$values), 1))
    ch <- chol(eig$vectors %*% (vals * t(eig$vectors)))
    flag$n <- flag$n + 1L
  }
  ch
}

# Log density of the rows of `Xc` (n x p) under N(mu, S) or, with `df`, the t
# distribution with location mu, scale matrix S and df degrees of freedom
.batchmixLogDensity <- function(Xc, mu, S, df, flag) {
  p <- length(mu)
  ch <- .batchmixCholSafe(S, flag)
  z <- backsolve(ch, t(Xc) - mu, transpose = TRUE)
  u <- colSums(z^2)
  half_logdet <- sum(log(diag(ch)))
  if (is.null(df)) {
    -0.5 * (p * log(2 * pi) + u) - half_logdet
  } else {
    lgamma((df + p) / 2) - lgamma(df / 2) - 0.5 * p * log(df * pi) - half_logdet -
      0.5 * (df + p) * log1p(u / df)
  }
}

.batchmixLogSumExp <- function(m) {
  mx <- apply(m, 1, max)
  mx + log(rowSums(exp(m - mx)))
}

.batchmixPointwiseChain <- function(ch, X, groups, fixed, index) {
  P <- ch$P
  K <- ch$K_max
  B <- ch$B
  N <- nrow(X)
  is_mvt <- identical(ch$type, "MVT")
  global <- is.null(ch$batch_weight_prior) || identical(ch$batch_weight_prior, "global")
  flag <- new.env()
  flag$n <- 0L

  ll <- matrix(0, length(index), N)
  for (i in seq_along(index)) {
    r <- index[i]
    mean_sum <- matrix(ch$mean_sum[, , r], nrow = P) # P x (K * B), cell k * B + b
    cov_comb <- matrix(ch$cov_comb[, , r], nrow = P) # P x (P * K * B)
    for (g in groups) {
      b <- g$batch
      o <- g$observed
      if (length(o) == 0) next
      items <- g$items
      w <- if (global) ch$weights[r, ] else ch$w_batch[b + 1, , r]
      comp <- matrix(0, length(items), K)
      for (k in seq_len(K)) {
        cell <- (k - 1) * B + b
        S <- cov_comb[, (cell * P + 1):(cell * P + P), drop = FALSE][o, o, drop = FALSE]
        df <- if (is_mvt) ch$t_df[r, k] else NULL
        comp[, k] <- log(w[k]) +
          .batchmixLogDensity(X[items, o, drop = FALSE], mean_sum[o, cell + 1], S, df, flag)
      }
      value <- .batchmixLogSumExp(comp)
      # an observed label is data: the joint density at that label
      lab <- fixed[items]
      if (any(lab)) {
        z <- ch$samples[r, items[lab]] + 1
        value[lab] <- comp[cbind(which(lab), z)]
      }
      ll[i, items] <- value
    }
  }
  list(ll = ll, n_regularised = flag$n)
}
