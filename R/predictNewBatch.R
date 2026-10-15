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
# Batch shift (m_b) and scale (S_b) each have a genuine partial-pooling/
# shrinkage prior (Gelman et al., 2013, BDA3, ch. 5) that can be estimated
# from the batches actually fitted, though - due to an identifiability
# constraint neither can escape - only ONE hyperparameter of each is ever
# free: m_b ~ N(0, delta_2 * lambda_2), with lambda_2 (= m_scale) ESTIMATED
# when sample_m_scale = TRUE; S_b ~ S_loc + InvGamma(rho, theta), with rho
# (the concentration; theta is pinned so the prior MEAN stays fixed at
# s_scale_prior_mean - see mvnSampler::sScaleConcentrationMetropolis()'s
# documentation for why the mean itself can never be freed) ESTIMATED when
# sample_s_scale = TRUE. When either flag was FALSE at training time, the
# corresponding predictive draw below falls back to the fixed, user-
# supplied hyperparameter exactly as before.

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
#
# Returns a list: `weight_draws` (n_saved x K, always present - the public
# result for the no-data case) plus, when weight_prior != "global" and
# there's more than one class, the auxiliary per-draw quantities
# predictNewBatch()'s X_new path needs to seed and then refine the same
# draw via composition sampling: `eta_known` (B x n_free x n_saved, the
# ORIGINAL batches' relabelling-safe ALR coordinates), `eta_star` (n_saved
# x n_free, the new batch's own seeded predictive draw - identical to what
# produced `weight_draws`), and the Rao-Blackwellised population
# hyperparameter(s) to hold fixed throughout prediction: `beta_hat`
# (n_free x n_saved, "gp") or `mu_hat`/`tau2_hat` (n_free x n_saved each,
# "partial_pooling").
.predictiveWeightDraws <- function(processed_chain, new_batch_coordinate = NULL) {
  K <- processed_chain$K_max
  B <- processed_chain$B
  weight_prior <- processed_chain$batch_weight_prior %||% "global"
  n_saved <- if (weight_prior == "global") {
    nrow(processed_chain$weights)
  } else {
    dim(processed_chain$w_batch)[3]
  }

  if (weight_prior == "gp" && is.null(new_batch_coordinate)) {
    stop("new_batch_coordinate is required when batch_weight_prior = 'gp'.")
  }
  if (weight_prior != "gp" && !is.null(new_batch_coordinate)) {
    warning("new_batch_coordinate is ignored unless batch_weight_prior = 'gp'.")
  }

  if (weight_prior == "global") {
    return(list(weight_draws = processed_chain$weights))
  }

  n_free <- K - 1
  weight_draws <- matrix(1, n_saved, K) # K == 1 edge case: trivially all mass on the one class
  if (n_free == 0) {
    return(list(weight_draws = weight_draws))
  }

  eta_star <- matrix(0, n_saved, n_free)
  eta_known_arr <- array(0, dim = c(B, n_free, n_saved))
  beta_hat <- mu_hat <- tau2_hat <- matrix(0, n_free, n_saved)
  mu_prior_sd <- processed_chain$pp_mu_prior_sd %||% 10.0

  for (t in seq_len(n_saved)) {
    wb <- processed_chain$w_batch[, , t]
    eta_known <- log(wb[, seq_len(n_free), drop = FALSE] / wb[, K]) # B x n_free, relabelling-safe
    eta_known_arr[, , t] <- eta_known

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
        beta_hat[j, t] <- beta_j
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
        mu_hat[j, t] <- post$mu
        tau2_hat[j, t] <- post$tau2
        eta_star[t, j] <- stats::rnorm(1, post$mu, sqrt(post$tau2))
      }
    }
  }

  D <- 1 + rowSums(exp(eta_star))
  weight_draws[, seq_len(n_free)] <- exp(eta_star) / D
  weight_draws[, K] <- 1 / D

  list(
    weight_draws = weight_draws, eta_known = eta_known_arr, eta_star = eta_star,
    beta_hat = beta_hat, mu_hat = mu_hat, tau2_hat = tau2_hat
  )
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

  # rho is always a per-iteration trace (constant unless sample_s_scale -
  # see batchSemiSupervisedMixtureModel.R); theta is kept in sync with it
  # (theta = s_scale_prior_mean * (rho - 1), fixing the prior MEAN of
  # S_b - S_loc - see sScaleConcentrationMetropolis()'s documentation for
  # why only the concentration, not the mean, is ever free), so it must be
  # recomputed per draw here rather than read as processed_chain$theta
  # (the ORIGINAL fixed value, only still correct when rho never moved).
  rho_draws <- processed_chain$rho
  theta_draws <- if (isTRUE(processed_chain$sample_s_scale)) {
    processed_chain$s_scale_prior_mean * (rho_draws - 1.0)
  } else {
    rep(processed_chain$theta, n_saved)
  }

  S_loc <- 1.0
  scale_draws <- S_loc + matrix(1 / stats::rgamma(P * n_saved,
    shape = rep(rho_draws, each = P), rate = rep(theta_draws, each = P)
  ), nrow = P, ncol = n_saved)

  list(shift = shift_draws, scale = scale_draws)
}

.predictNewBatchDriverName <- function(type) {
  switch(type,
    MVN = "predictNewBatchMVN",
    MVT = "predictNewBatchMVT",
    MVN_LKJ = "predictNewBatchMVNSeparationStrategy",
    MVN_MIXED = "predictNewBatchMVNMixed",
    stop("Type not recognised. Please use one of 'MVN', 'MVT', 'MVN_LKJ' or 'MVN_MIXED'.")
  )
}

#' @title Predict a new batch
#' @description Posterior-predictive quantities for a batch not seen
#' during fitting, for any model \code{type} and any \code{batch_weight_prior}.
#' With no new data (\code{X_new = NULL}), draws the new batch's own class
#' weights and location/scale shift straight from their fitted priors -
#' for \code{batch_weight_prior = "gp"} this is a genuine Gaussian process
#' prediction at \code{new_batch_coordinate} (e.g. a future collection
#' time). With \code{X_new} supplied, additionally classifies that data by
#' composition sampling: for each retained posterior draw, the new batch's
#' own shift/scale/weight and item allocations are drawn from their
#' conditional posterior GIVEN that draw's (fixed) cluster parameters and
#' \code{X_new}'s own likelihood, by literally reusing the training
#' sampler's own Metropolis/Gibbs updates restricted to the one new batch
#' (\code{sampler::predict_mode}, src/sampler.h) - see the file-level
#' comment in \code{R/predictNewBatch.R} for the full statistical
#' justification (composition sampling for a posterior predictive
#' distribution, Rubin 1987; BDA3 sec. 1.10). Batch shift/scale each draw
#' from a genuine cross-batch partial-pooling prior when the original fit
#' set \code{sample_m_scale}/\code{sample_s_scale = TRUE} respectively (a
#' fixed, non-adaptive prior otherwise - see that same comment).
#' @param processed_chain Output of \code{\link{processMCMCChain}} (a
#' single chain, burn-in applied and relabelled).
#' @param X The original training data passed to the function that
#' produced \code{processed_chain} (same convention as
#' \code{\link{continueChain}}'s \code{X} argument).
#' @param batch_vec The ORIGINAL per-item batch labels for \code{X} (same
#' convention as \code{\link{continueChain}}'s \code{batch_vec} argument -
#' the fit object does not store this itself). Required, and only used,
#' when \code{X_new} is supplied.
#' @param new_batch_coordinate The new batch's 1-D coordinate (e.g. its
#' collection time), required when \code{processed_chain$batch_weight_prior == "gp"}
#' and ignored otherwise.
#' @param X_new The new batch's own data (items in rows, same \code{P}
#' columns as \code{X}), or \code{NULL} (the default) for the no-data
#' prior/GP-predictive case.
#' @param fixed_new Binary vector, one entry per row of \code{X_new}: 1 for
#' any item whose label is already known (semi-supervised), 0 otherwise.
#' Defaults to all-unknown. Ignored if \code{X_new} is \code{NULL}.
#' @param labels_new_init Known/initial labels for \code{X_new}'s items
#' (1-indexed, matching this package's usual convention); only entries
#' where \code{fixed_new == 1} matter. Defaults to all-1 (arbitrary,
#' since every unfixed item's label is redrawn from scratch regardless of
#' its starting value). Ignored if \code{X_new} is \code{NULL}.
#' @param n_draws Number of retained posterior iterations from
#' \code{processed_chain} to use for composition sampling; \code{NULL}
#' (the default) uses every retained iteration. Lower this to trade
#' predictive-sample size for speed on a large chain. Ignored if
#' \code{X_new} is \code{NULL}.
#' @param n_pred_iter Number of Metropolis/Gibbs sweeps run per posterior
#' draw. Ignored if \code{X_new} is \code{NULL}.
#' @param pred_burn Sweeps discarded (per draw) before recording starts;
#' defaults to \code{floor(n_pred_iter / 2)}. Ignored if \code{X_new} is
#' \code{NULL}.
#' @param pred_thin Thinning factor applied to the recorded sweeps.
#' Ignored if \code{X_new} is \code{NULL}.
#' @param column_type,censor_code Only for \code{processed_chain$type == "MVN_MIXED"}:
#' the SAME \code{column_type}/\code{censor_code} originally passed to
#' \code{\link{batchSemiSupervisedMixtureModel}} (again not stored on the
#' fit object - same convention as \code{X}/\code{batch_vec}). Required,
#' and only used, when \code{X_new} is supplied and the type is
#' \code{"MVN_MIXED"}.
#' @param censor_code_new Only for \code{type == "MVN_MIXED"}: \code{X_new}'s
#' own censoring indicator, same convention as \code{censor_code}. Defaults
#' to no censoring (all-zero) if \code{X_new} has any censored/binary
#' columns but this is left \code{NULL}.
#' @return A named list. Always present:
#' \itemize{
#'   \item \code{weight_draws}: \code{n_saved x K} matrix, one predictive
#'   class-weight draw per retained posterior iteration (no-data
#'   prior/GP-predictive draw; also used to seed the composition-sampling
#'   draw when \code{X_new} is supplied).
#'   \item \code{shift_draws}, \code{scale_draws}: \code{P x n_saved}
#'   matrices, the new batch's predictive location shift/scale.
#'   \item \code{weight_est}, \code{shift_est}, \code{scale_est}: median
#'   point estimates.
#'   \item \code{mean_sum_est}: \code{P x K} matrix, \code{mean_est + shift_est}.
#' }
#' Additionally, when \code{X_new} is supplied: \code{label_new} (a matrix,
#' one row per composition-sampling draw, one column per row of
#' \code{X_new}), \code{allocation_probability} (\code{N_new x K} point
#' estimate, as \code{\link{calcAllocProb}}), \code{prob}/\code{pred} (the
#' usual per-item probability/class point estimates), and
#' \code{shift_new_draws}/\code{scale_new_draws}/\code{weight_new_draws}
#' (the composition-sampling-refined versions of \code{shift_draws}/
#' \code{scale_draws}/\code{weight_draws}, now conditioned on \code{X_new}
#' as well as the fitted model).
#' @export
predictNewBatch <- function(processed_chain, X, batch_vec = NULL, new_batch_coordinate = NULL, X_new = NULL,
                             fixed_new = NULL, labels_new_init = NULL,
                             n_draws = NULL, n_pred_iter = 200, pred_burn = NULL, pred_thin = 1,
                             column_type = NULL, censor_code = NULL, censor_code_new = NULL) {
  if (is.null(processed_chain$mean_est)) {
    stop("processed_chain must be the output of processMCMCChain() (burn-in applied, point estimates computed).")
  }
  if (!is.null(X_new)) {
    if (isTRUE(processed_chain$include_interaction)) {
      stop(
        "predictNewBatch() with X_new does not support include_interaction = TRUE fits: ",
        "the fitted batch x cluster interaction term would be silently dropped for every ",
        "batch (not just the new one), not merely omitted for the new batch alone."
      )
    }
    if (is.null(batch_vec)) {
      stop("batch_vec (the original per-item batch labels) is required when X_new is supplied.")
    }
  }

  weights_info <- .predictiveWeightDraws(processed_chain, new_batch_coordinate)
  weight_draws <- weights_info$weight_draws
  shift_scale <- .predictiveShiftScaleDraws(processed_chain, X)

  weight_est <- apply(weight_draws, 2, stats::median)
  shift_est <- apply(shift_scale$shift, 1, stats::median)
  scale_est <- apply(shift_scale$scale, 1, stats::median)

  out <- list(
    weight_draws = weight_draws,
    shift_draws = shift_scale$shift,
    scale_draws = shift_scale$scale,
    weight_est = weight_est,
    shift_est = shift_est,
    scale_est = scale_est,
    mean_sum_est = processed_chain$mean_est + shift_est
  )

  if (is.null(X_new)) {
    return(out)
  }

  N_new <- nrow(X_new)
  fixed_new <- fixed_new %||% rep(0L, N_new)
  labels_new_init <- labels_new_init %||% rep(1L, N_new)

  n_saved_total <- dim(processed_chain$means)[3]
  draw_idx <- if (is.null(n_draws)) {
    seq_len(n_saved_total)
  } else {
    round(seq(1, n_saved_total, length.out = min(n_draws, n_saved_total)))
  }
  if (is.null(pred_burn)) pred_burn <- floor(n_pred_iter / 2)

  type <- processed_chain$type
  K <- processed_chain$K_max
  B <- processed_chain$B
  weight_prior <- processed_chain$batch_weight_prior %||% "global"
  weight_prior_type <- switch(weight_prior, global = 0L, partial_pooling = 1L, gp = 2L)
  n_free <- max(K - 1, 0)

  batch_coordinates_new <- if (weight_prior == "gp") {
    c(processed_chain$batch_coordinates, new_batch_coordinate)
  } else {
    numeric(B + 1)
  }

  eta_alr_init_draws <- array(0, dim = c(B + 1, n_free, length(draw_idx)))
  if (weight_prior != "global" && n_free > 0) {
    for (i in seq_along(draw_idx)) {
      t <- draw_idx[i]
      # matrix(..., ncol = n_free), not a bare [, , t] slice: with n_free == 1
      # (K == 2) that slice silently drops to a plain B-length vector,
      # breaking rbind()'s column alignment below.
      eta_known_t <- matrix(weights_info$eta_known[, , t], ncol = n_free)
      eta_star_t <- matrix(weights_info$eta_star[t, ], nrow = 1, ncol = n_free)
      eta_alr_init_draws[, , i] <- rbind(eta_known_t, eta_star_t)
    }
  }

  batch_vec_0idx <- as.integer(as.factor(batch_vec)) - 1L
  if (max(batch_vec_0idx) + 1L != B) {
    stop("batch_vec does not have exactly B = ", B, " distinct values matching processed_chain.")
  }

  args <- list(
    X = X, X_new = X_new, K = K, B = B,
    batch_vec = batch_vec_0idx,
    fixed_new = as.integer(fixed_new),
    labels_new_init = as.integer(labels_new_init) - 1L,
    label_draws = t(processed_chain$samples[draw_idx, , drop = FALSE]),
    means_draws = processed_chain$means[, , draw_idx, drop = FALSE],
    cov_draws = processed_chain$covariance[, , draw_idx, drop = FALSE],
    batch_shift_draws = processed_chain$batch_shift[, , draw_idx, drop = FALSE],
    batch_scale_draws = processed_chain$batch_scale[, , draw_idx, drop = FALSE],
    shift_new_init = shift_scale$shift[, draw_idx, drop = FALSE],
    scale_new_init = shift_scale$scale[, draw_idx, drop = FALSE],
    m_scale = processed_chain$m_scale %||% 0.01,
    lambda_2_draws = if (isTRUE(processed_chain$sample_m_scale)) processed_chain$lambda_2[draw_idx] else rep(0, length(draw_idx)),
    sample_m_scale = isTRUE(processed_chain$sample_m_scale),
    # rho is always a per-iteration trace (see processMCMCChain.R's note);
    # its first retained value is only ever used as the constructor's
    # initial value here - every draw's own rho_draws(t) overwrites it
    # before matrixCombinations() runs, whether or not sample_s_scale was
    # actually TRUE at training time (see the 4 predictNewBatch*.cpp
    # drivers' `if (sample_s_scale) { ... }` guard).
    rho = processed_chain$rho[1],
    theta = processed_chain$theta,
    rho_draws = processed_chain$rho[draw_idx],
    s_scale_prior_mean = processed_chain$s_scale_prior_mean %||% 1.0,
    sample_s_scale = isTRUE(processed_chain$sample_s_scale),
    weight_prior_type = weight_prior_type,
    weights_draws = if (weight_prior == "global") weight_draws[draw_idx, , drop = FALSE] else matrix(0, length(draw_idx), K),
    eta_alr_init_draws = eta_alr_init_draws,
    gp_beta_draws = if (weight_prior == "gp") weights_info$beta_hat[, draw_idx, drop = FALSE] else matrix(0, n_free, length(draw_idx)),
    pp_mu_draws = if (weight_prior == "partial_pooling") weights_info$mu_hat[, draw_idx, drop = FALSE] else matrix(0, n_free, length(draw_idx)),
    pp_tau2_draws = if (weight_prior == "partial_pooling") weights_info$tau2_hat[, draw_idx, drop = FALSE] else matrix(0, n_free, length(draw_idx)),
    gp_tau2_draws = if (weight_prior == "gp") processed_chain$gp_tau2[draw_idx] else rep(1, length(draw_idx)),
    gp_length_scale_draws = if (weight_prior == "gp") processed_chain$gp_length_scale[draw_idx] else rep(1, length(draw_idx)),
    batch_coordinates_new = batch_coordinates_new,
    eta_proposal_window = processed_chain$final_eta_proposal_window %||% processed_chain$eta_proposal_window %||% 0.1,
    m_proposal_window = processed_chain$final_m_proposal_window %||% processed_chain$m_proposal_window %||% 0.5**2,
    S_proposal_window = if (!is.null(processed_chain$final_S_proposal_window)) {
      1.0 / processed_chain$final_S_proposal_window
    } else {
      processed_chain$S_proposal_window %||% 0.01
    },
    n_pred_iter = n_pred_iter, pred_burn = pred_burn, pred_thin = pred_thin
  )

  if (type == "MVT") {
    args$t_df_draws <- processed_chain$t_df[draw_idx, , drop = FALSE] # n_draws x K, matching predictNewBatchMVT()'s t_df_draws.row(t)
    args$t_df_proposal_window <- processed_chain$t_df_proposal_window %||% 0.015
    # matches predictNewBatchMVT()'s parameter order (t_df_draws sits right
    # after cov_draws; t_df_proposal_window right after S_proposal_window) -
    # do.call() matches by name, so list order here doesn't actually matter,
    # but keeping this comment next to both insertions for anyone reading it.
  }
  if (type == "MVN_MIXED") {
    if (is.null(column_type) || is.null(censor_code)) {
      stop("column_type and censor_code (the same ones passed to the original fit) are required for type = 'MVN_MIXED'.")
    }
    args$column_type <- as.integer(column_type)
    args$censor_code <- matrix(as.integer(censor_code), nrow(X), ncol(X))
    args$censor_code_new <- if (is.null(censor_code_new)) {
      matrix(0L, N_new, ncol(X_new))
    } else {
      matrix(as.integer(censor_code_new), N_new, ncol(X_new))
    }
  }

  raw <- do.call(get(.predictNewBatchDriverName(type)), args)

  # raw$label_new/alloc_new are 0-indexed (the C++ layer's convention);
  # +1 to match this package's usual 1-indexed class labels everywhere else
  # (e.g. calcAllocProb()/predictClass()'s output).
  label_new <- raw$label_new + 1L
  alloc_prob <- apply(raw$alloc_new, c(1, 2), median)

  out$label_new <- label_new
  out$allocation_probability <- alloc_prob
  out$prob <- apply(alloc_prob, 1, max)
  out$pred <- apply(alloc_prob, 1, which.max)
  out$shift_new_draws <- raw$shift_new
  out$scale_new_draws <- raw$scale_new
  out$weight_new_draws <- t(raw$weight_new)
  out$shift_est <- apply(raw$shift_new, 1, stats::median)
  out$scale_est <- apply(raw$scale_new, 1, stats::median)
  out$weight_est <- apply(raw$weight_new, 1, stats::median)
  out$mean_sum_est <- processed_chain$mean_est + out$shift_est

  out
}
