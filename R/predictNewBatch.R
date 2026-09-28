#!/usr/bin/Rscript
# Predicting a new batch (see ?predictNewBatch).
#
# Two genuinely different predictive problems share this name:
#
#  (i) No data yet for the new batch (e.g. "what will next month's batch
#      look like"): draw the new batch's own weight/shift/scale straight
#      from their fitted priors - for `batch_weight_prior = "gp"`, this is
#      the standard Gaussian process conditional-predictive at a new
#      coordinate (Rasmussen & Williams, 2006, "Gaussian Processes for
#      Machine Learning," MIT Press, Eq. 2.19); for "partial_pooling", a
#      draw from the estimated population distribution of batch weights
#      (the same "predict a new group" logic as Gelman & Hill, 2007,
#      "Data Analysis Using Regression and Multilevel/Hierarchical
#      Models," ch. 12, and McElreath, 2020, "Statistical Rethinking," 2nd
#      ed., ch. 13-14); for "global" the shared weight is simply reused.
#      Batch shift/scale have no such cross-batch hierarchical mean under
#      any weight-prior choice (see the block comment below), so their
#      predictive draw is always just the fitted prior.
#
#  (ii) New (unlabelled) data has arrived under the new batch: its own
#      shift/scale/weight are inferred from that data while the already-
#      fitted cluster parameters stay fixed, then its items are
#      classified - a composition-sampling posterior predictive (Rubin,
#      1987, "Multiple Imputation for Nonresponse in Surveys," Wiley,
#      sec. 2.3; Gelman, Carlin, Stern, Dunson, Vehtari & Rubin, 2013,
#      "Bayesian Data Analysis," 3rd ed. (BDA3), sec. 1.10, for drawing
#      from the posterior predictive distribution by composition):
#      simulate from p(new-batch parameters | posterior draw s, new data)
#      once per retained posterior draw s, via the same conditional
#      updates the training sampler itself uses for one batch (see
#      `sampler::predict_mode` in src/sampler.h), so the classification
#      accounts for both the new batch's own likelihood AND the fitted
#      model's remaining posterior uncertainty.
#
# Batch shift (m_b) and scale (S_b) are NOT pooled across batches the way
# weights are: m_b ~ N(0, delta_2 * lambda_2) with lambda_2 (= m_scale)
# ESTIMATED from all batches jointly when sample_m_scale = TRUE (a
# genuine partial-pooling/shrinkage prior, Gelman et al., 2013, BDA3,
# ch. 5), but S_b ~ S_loc + InvGamma(rho, theta) with rho/theta fixed,
# user-supplied hyperparameters that are never re-estimated from the
# batches actually fitted. This asymmetry means a new batch's shift
# predictive correctly reflects what the fitted batches' own shifts
# looked like, but its scale predictive currently falls back to a
# generic, non-adaptive prior regardless of the fitted batches' own
# scales - a real limitation, not addressed here (it would require adding
# a new hierarchical hyperparameter to the core training model, well
# beyond the scope of a predictive function built on top of the existing
# fit) but worth knowing about when interpreting a new batch's predicted
# scale.

# Local null-coalescing helper (see R/continueChain.R for the identical,
# independently-defined copy - this file does not depend on that one).
`%||%` <- function(a, b) if (is.null(a)) b else a

# Mean-impute then take the mean diagonal (co)variance, exactly mirroring
# mvnSampler's C++ constructor (imputeColumnMeans() + arma::cov(), which
# uses the same N-1 sample-covariance normalisation as R's cov()) - see
# genericFunctions.h. This is only ever used as a weakly-informative
# empirical-Bayes SCALE for the shift prior (delta_2), not a
# posterior-inferred quantity, so this one-off recomputation from the
# original data is exact, not an approximation.
.batchShiftDelta2 <- function(X) {
  X_imp <- X
  for (p in seq_len(ncol(X))) {
    col <- X_imp[, p]
    finite <- is.finite(col)
    if (!all(finite)) col[!finite] <- mean(col[finite])
    X_imp[, p] <- col
  }
  mean(diag(stats::cov(X_imp)))
}

# Exact GLS conjugate posterior mean of the GP intercept beta_j given one
# free ALR coordinate's values `eta_j` (B-vector) and the GP covariance's
# lower-triangular Cholesky factor `gp_chol` (Sigma = gp_chol %*% t(gp_chol))
# - the same formula as sampler::updateGPWeights()'s Gibbs step in
# src/sampler.cpp, evaluated here as a posterior MEAN rather than a fresh
# Gibbs draw (Rao-Blackwellisation: an exact, lower-variance, unbiased
# summary of the same conditional distribution - Casella & Robert, 1996,
# "Rao-Blackwellisation of sampling schemes," Biometrika 83(1)). Used
# instead of the raw sampled gp_beta because gp_beta is indexed by ALR
# *coordinate*, not by class, and so cannot be safely reused once
# relabelChain() has permuted classes (see the "relabelling safety" note
# in .predictiveWeightDraws() below) - recomputing it from the
# already-relabelled eta is exact, not an approximation, since this
# formula's only inputs are eta_j and gp_chol (which does not depend on
# class labels at all).
.gpBetaPosteriorMean <- function(eta_j, gp_chol, mu_prior_sd) {
  ones_B <- rep(1, length(eta_j))
  L_inv_ones <- forwardsolve(gp_chol, ones_B)
  L_inv_eta <- forwardsolve(gp_chol, eta_j)
  precision_beta <- sum(L_inv_ones^2) + 1 / mu_prior_sd^2
  sum(L_inv_ones * L_inv_eta) / precision_beta
}

# Self-consistent (mu_j, tau2_j) for the partial-pooling prior, given one
# free ALR coordinate's values `eta_j`. Unlike the GP intercept above,
# mu_j and tau2_j are each other's exact conjugate posterior MEAN only
# conditional on the other (sampler::updatePartialPoolingWeights()'s two
# separate Gibbs steps, not a single closed-form joint), so this iterates
# the two conditional-mean updates to a fixed point (a handful of
# iterations; B is always small) rather than picking one arbitrarily -
# the semi-conjugate analogue of the exact GP case above (the prior
# mu_j ~ N(0, mu_prior_sd^2) does not scale with tau2_j, so the joint
# posterior is not itself a standard Normal-InverseGamma).
.ppMuTau2PosteriorMean <- function(eta_j, mu_prior_sd, tau2_prior_shape, tau2_prior_rate,
                                    n_fixed_point_iter = 10) {
  B <- length(eta_j)
  mu_j <- mean(eta_j)
  tau2_j <- tau2_prior_rate / max(tau2_prior_shape - 1, 1e-6)
  for (iter in seq_len(n_fixed_point_iter)) {
    var_post_mu <- 1 / (1 / mu_prior_sd^2 + B / tau2_j)
    mu_j <- var_post_mu * sum(eta_j) / tau2_j
    shape_post <- tau2_prior_shape + B / 2
    rate_post <- tau2_prior_rate + 0.5 * sum((eta_j - mu_j)^2)
    tau2_j <- rate_post / (shape_post - 1)
  }
  list(mu = mu_j, tau2 = tau2_j)
}

# One posterior-predictive draw of a new batch's ALR coordinates (and
# hence its class weight vector), for every retained iteration of
# `processed_chain` at once. NULL `new_batch_coordinate` is only valid
# for "global"/"partial_pooling" (batches have no coordinate to place a
# new one relative to); "gp" requires it.
#
# Relabelling safety: relabelChain() (run inside processMCMCChain())
# correctly permutes `w_batch`'s class (column) order every iteration,
# but does NOT (and, given the ALR transform's dependence on a fixed
# pivot class, cannot cheaply) permute the separately-stored `eta_alr`/
# `gp_beta`/`pp_mu`/`pp_tau2` traces, which are indexed by ALR
# *coordinate* under the model's ORIGINAL, possibly label-switched class
# ordering. Using those fields directly here would silently mix labels
# whenever the chain actually switched. This function therefore never
# reads them: every ALR quantity it needs (`eta_j`, and hence `beta_j`/
# `mu_j`/`tau2_j`) is recomputed from the already-correctly-relabelled
# `w_batch` instead (see .gpBetaPosteriorMean()/.ppMuTau2PosteriorMean()
# above) - exact for beta_j, and a principled fixed-point summary for
# (mu_j, tau2_j). gp_tau2/gp_length_scale are single values shared across
# every coordinate (not class-indexed), so are unaffected by relabelling
# and are used as sampled.
.predictiveWeightDraws <- function(processed_chain, new_batch_coordinate = NULL) {
  K <- processed_chain$K_max
  B <- processed_chain$B
  weight_prior <- processed_chain$batch_weight_prior %||% "global"
  n_saved <- if (weight_prior == "global") {
    nrow(processed_chain$weights)
  } else {
    dim(processed_chain$w_batch)[3]
  }

  if (weight_prior == "global") {
    return(processed_chain$weights)
  }

  if (weight_prior == "gp" && is.null(new_batch_coordinate)) {
    stop("new_batch_coordinate is required when batch_weight_prior = 'gp'.")
  }
  if (weight_prior != "gp" && !is.null(new_batch_coordinate)) {
    warning("new_batch_coordinate is ignored unless batch_weight_prior = 'gp'.")
  }

  n_free <- K - 1
  weight_draws <- matrix(1, n_saved, K) # K == 1 edge case: trivially all mass on the one class
  if (n_free == 0) {
    return(weight_draws)
  }

  eta_star <- matrix(0, n_saved, n_free)
  mu_prior_sd <- processed_chain$pp_mu_prior_sd %||% 10.0

  for (t in seq_len(n_saved)) {
    wb <- processed_chain$w_batch[, , t]
    eta_known <- log(wb[, seq_len(n_free), drop = FALSE] / wb[, K]) # B x n_free, relabelling-safe

    if (weight_prior == "gp") {
      gp_tau2 <- processed_chain$gp_tau2[t]
      gp_length_scale <- processed_chain$gp_length_scale[t]
      coords <- c(processed_chain$batch_coordinates, new_batch_coordinate)
      full_cov <- maternKernel32(coords, gp_tau2, gp_length_scale, jitter = 1e-8)
      Sigma_BB <- full_cov[seq_len(B), seq_len(B), drop = FALSE]
      k_star <- full_cov[B + 1, seq_len(B)]
      k_star_star <- full_cov[B + 1, B + 1]
      L <- t(chol(Sigma_BB))

      for (j in seq_len(n_free)) {
        e <- eta_known[, j]
        beta_j <- .gpBetaPosteriorMean(e, L, mu_prior_sd)
        solved <- forwardsolve(L, e - beta_j)
        solved <- backsolve(t(L), solved)
        mean_star <- beta_j + sum(k_star * solved)
        Sigma_BB_inv_kstar <- backsolve(t(L), forwardsolve(L, k_star))
        var_star <- max(k_star_star - sum(k_star * Sigma_BB_inv_kstar), 1e-8)
        eta_star[t, j] <- stats::rnorm(1, mean_star, sqrt(var_star))
      }
    } else { # partial_pooling
      tau2_prior_shape <- processed_chain$pp_tau2_prior_shape %||% 2.0
      tau2_prior_rate <- processed_chain$pp_tau2_prior_rate %||% 1.0
      for (j in seq_len(n_free)) {
        e <- eta_known[, j]
        post <- .ppMuTau2PosteriorMean(e, mu_prior_sd, tau2_prior_shape, tau2_prior_rate)
        eta_star[t, j] <- stats::rnorm(1, post$mu, sqrt(post$tau2))
      }
    }
  }

  D <- 1 + rowSums(exp(eta_star))
  weight_draws[, seq_len(n_free)] <- exp(eta_star) / D
  weight_draws[, K] <- 1 / D
  weight_draws
}

# Posterior-predictive draws of a new batch's own shift (m*) and scale
# (S*), one draw per retained iteration of `processed_chain`, straight
# from their fitted priors (see the block comment above for why these
# have no cross-batch pooled mean to draw towards, unlike the weight
# side). `X` is the ORIGINAL training data (same convention as
# continueChain()'s `X` argument): the fit object does not store it, and
# it is needed here only to recompute delta_2, the empirical-Bayes shift
# prior scale (see .batchShiftDelta2()).
.predictiveShiftScaleDraws <- function(processed_chain, X) {
  P <- processed_chain$P
  n_saved <- dim(processed_chain$means)[3]

  delta_2 <- .batchShiftDelta2(X)
  sample_m_scale <- isTRUE(processed_chain$sample_m_scale)
  lambda_2_draws <- if (sample_m_scale) {
    processed_chain$lambda_2
  } else {
    rep(processed_chain$m_scale, n_saved)
  }

  shift_sd <- sqrt(delta_2 * lambda_2_draws)
  shift_draws <- matrix(stats::rnorm(P * n_saved, mean = 0, sd = rep(shift_sd, each = P)),
    nrow = P, ncol = n_saved
  )

  S_loc <- 1.0
  scale_draws <- S_loc + matrix(1 / stats::rgamma(P * n_saved,
    shape = processed_chain$rho, rate = processed_chain$theta
  ), nrow = P, ncol = n_saved)

  list(shift = shift_draws, scale = scale_draws)
}

#' @title Predict a new batch
#' @description Posterior-predictive quantities for a batch not seen
#' during fitting, for any model \code{type} and any \code{batch_weight_prior}.
#' With no new data (\code{X_new = NULL}), draws the new batch's own class
#' weights, and location/scale shift, straight from their fitted priors -
#' for \code{batch_weight_prior = "gp"} this is a genuine Gaussian process
#' prediction at \code{new_batch_coordinate} (e.g. a future collection
#' time); see the file-level comment in \code{R/predictNewBatch.R} for the
#' full statistical justification, including why \code{"gp"}/\code{"partial_pooling"}
#' weights have a real cross-batch predictive distribution but batch scale
#' currently does not (a known asymmetry in the underlying model, not
#' something this function can fix). \code{X_new} (classifying actual new
#' data arriving under the new batch) is not yet supported by this
#' function.
#' @param processed_chain Output of \code{\link{processMCMCChain}} (a
#' single chain, burn-in applied and relabelled).
#' @param X The original training data passed to the function that
#' produced \code{processed_chain} (same convention as
#' \code{\link{continueChain}}'s \code{X} argument) - needed only to
#' recompute the batch-shift prior's empirical-Bayes scale.
#' @param new_batch_coordinate The new batch's 1-D coordinate (e.g. its
#' collection time), required when \code{processed_chain$batch_weight_prior == "gp"}
#' and ignored otherwise.
#' @param X_new Not yet implemented; must be \code{NULL}.
#' @return A named list:
#' \itemize{
#'   \item \code{weight_draws}: \code{n_saved x K} matrix, one predictive
#'   class-weight draw per retained posterior iteration.
#'   \item \code{shift_draws}, \code{scale_draws}: \code{P x n_saved}
#'   matrices, the new batch's predictive location shift/scale.
#'   \item \code{weight_est}, \code{shift_est}, \code{scale_est}: median
#'   point estimates (vectors of length \code{K}/\code{P}/\code{P}).
#'   \item \code{mean_sum_est}: \code{P x K} matrix, \code{mean_est + shift_est}
#'   - the new batch's predicted, batch-shifted cluster locations.
#' }
#' @export
predictNewBatch <- function(processed_chain, X, new_batch_coordinate = NULL, X_new = NULL) {
  if (!is.null(X_new)) {
    stop(
      "Predicting a new batch's own data (X_new) is not yet implemented; ",
      "only the no-data (prior/GP-predictive) case is currently supported."
    )
  }
  if (is.null(processed_chain$mean_est)) {
    stop("processed_chain must be the output of processMCMCChain() (burn-in applied, point estimates computed).")
  }

  weight_draws <- .predictiveWeightDraws(processed_chain, new_batch_coordinate)
  shift_scale <- .predictiveShiftScaleDraws(processed_chain, X)

  weight_est <- apply(weight_draws, 2, stats::median)
  shift_est <- apply(shift_scale$shift, 1, stats::median)
  scale_est <- apply(shift_scale$scale, 1, stats::median)

  list(
    weight_draws = weight_draws,
    shift_draws = shift_scale$shift,
    scale_draws = shift_scale$scale,
    weight_est = weight_est,
    shift_est = shift_est,
    scale_est = scale_est,
    mean_sum_est = processed_chain$mean_est + shift_est
  )
}
