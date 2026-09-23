#' @title Simulate from the prior predictive distribution
#' @description Draws parameters from the model's priors (never conditioning
#' on any observed cluster/batch likelihood) and simulates full replicate
#' datasets from them, for a prior predictive check - "does the model, before
#' seeing any data, generate datasets that look like the kind of thing we
#' expect to see?" (Gelman, Carlin, Stern, Dunson, Vehtari & Rubin, 2013,
#' \emph{Bayesian Data Analysis}, 3rd ed., Section 6.1 and 6.3-6.4; Gelman et
#' al., 2020, "Bayesian Workflow", arXiv:2011.01808, Section 2.1). This is
#' the standard first check in a Bayesian workflow, done before the model
#' ever touches the real clustering/batch structure: if simulated datasets
#' routinely look absurd (implausible scales, degenerate clusters, batch
#' shifts far larger than any plausible batch effect), that is a sign the
#' prior - not yet the data or the likelihood - needs reconsidering.
#'
#' The priors used here are empirical-Bayes (their location/scale are
#' derived from \code{X}'s own summary statistics, exactly as the real
#' sampler's constructor computes them - see \code{src/mvnSampler.cpp} /
#' \code{src/mvnSamplerSeparationStrategy.cpp}), so this is not a check of
#' a data-independent prior; it checks whether the empirical-Bayes
#' construction itself produces sensible simulated data at this \code{K}/
#' \code{B}, before any label or cluster/batch structure has been fit.
#'
#' \strong{'MVN', 'MVT', 'MVN_LKJ' and 'MVN_MIXED' are supported.} For
#' \code{type = "MVN_MIXED"}, every column simulates the same underlying
#' latent Gaussian draw as \code{"MVN_LKJ"} (the model \code{"MVN_MIXED"}
#' shares with it), then binary/probit-linked columns
#' (\code{column_type == 1}) are thresholded at 0 (Albert & Chib, 1993) to
#' produce the replicated \{0, 1\} observation - matching \code{sigma}/
#' \code{S} being fixed at 1 for those columns in \code{sampleCovPrior()}/
#' \code{sampleSPrior()} (\code{src/mvnSamplerMixed.cpp}), the standard
#' probit-identification device. Every replicate is fully observed (no
#' \code{NA}/censoring), including at cells missing or censored in the real
#' data: the model's own assumption is that those cells are draws from
#' exactly the same distribution as every other cell (ignorable
#' missingness; censoring only affects what is *recorded*, not the
#' underlying generative draw), so simulating them like any other cell is
#' the model-consistent replicate. \code{\link{plotPredictiveCheck}}'s
#' \code{censor_code} argument instead excludes the real data's censored
#' cells from the *comparison*, since a recorded censoring bound is not the
#' true value to compare a free replicate draw against - see its
#' documentation.
#'
#' @param X Data matrix (items in rows), used only to derive the
#' empirical-Bayes prior location/scale (its own values are never reused
#' directly) - the same summary statistics
#' (\code{colMeans(X)}/\code{cov(X)}) that \code{batchSemiSupervisedMixtureModel()}
#' computes internally when \code{type} is 'MVN', 'MVT', 'MVN_LKJ' or
#' 'MVN_MIXED'. Missing entries (\code{NA}/\code{NaN}) are allowed here too -
#' they are replaced by their column mean only for this hyperparameter
#' derivation (mirroring the C++ samplers' own empirical-Bayes setup
#' exactly), never imputed via the real per-sweep data augmentation used in
#' an actual fit, since a prior predictive draw does not condition on any
#' observed data at all. For \code{type = "MVN_MIXED"}, binary columns'
#' raw \{0, 1\} values are used as-is for this derivation too (again
#' mirroring the C++ constructor exactly), not any latent-scale transform.
#' @param batch_vec Observed batch label for each row of X.
#' @param K The number of clusters to simulate.
#' @param type One of 'MVN', 'MVT', 'MVN_LKJ', 'MVN_MIXED' (see
#' \code{\link{batchSemiSupervisedMixtureModel}} for what each means).
#' @param column_type Only used if \code{type = "MVN_MIXED"}: a P-vector,
#' \code{0} for a continuous column, \code{1} for a binary/probit column -
#' see \code{\link{batchSemiSupervisedMixtureModel}}. Required for that type.
#' @param alpha Symmetric Dirichlet concentration for the prior cluster
#' weights (ignored if \code{concentration} is given).
#' @param concentration K-vector of Dirichlet concentrations for the prior
#' cluster weights; defaults to \code{rep(alpha, K)}.
#' @param m_scale The (fixed) batch-shift prior scale hyperparameter - see
#' \code{batchSemiSupervisedMixtureModel(sample_m_scale = ...)}; this
#' function always treats it as fixed, i.e. as if \code{sample_m_scale =
#' FALSE} (the extra InvGamma hyperprior on \code{m_scale} itself is not
#' simulated).
#' @param rho,theta Shape/rate of the batch-scale prior.
#' @param eta LKJ concentration parameter; only used if \code{type} is
#' 'MVN_LKJ' or 'MVN_MIXED'.
#' @param t_df_shape,t_df_rate,t_df_loc Hyperparameters of the (shifted
#' Gamma) prior on the cluster degrees of freedom; only used if \code{type}
#' is 'MVT'. Defaults match \code{mvtSampler}'s own defaults
#' (\code{psi = 2, chi = 0.1, t_loc = 2}).
#' @param n_datasets Number of independent prior predictive datasets to
#' draw.
#' @return A list of length \code{n_datasets}; each element is a list with
#' \code{X} (the simulated, fully-observed N x P data matrix), \code{labels}
#' (the simulated N-vector of prior cluster draws, 0-indexed) and
#' \code{params} (the drawn mu/cov/m/S/weights, and t_df if \code{type ==
#' "MVT"}).
#' @seealso \code{\link{simulatePosteriorPredictive}},
#' \code{\link{plotPredictiveCheck}}
#' @export
simulatePriorPredictive <- function(X,
                                    batch_vec,
                                    K,
                                    type = c("MVN", "MVT", "MVN_LKJ", "MVN_MIXED"),
                                    column_type = NULL,
                                    alpha = 1,
                                    concentration = NULL,
                                    m_scale = 0.01,
                                    rho = 3.0,
                                    theta = 1.0,
                                    eta = 1.0,
                                    t_df_shape = 2.0,
                                    t_df_rate = 0.1,
                                    t_df_loc = 2.0,
                                    n_datasets = 1) {
  type <- match.arg(type)
  is_mixed <- type == "MVN_MIXED"

  if (!is.matrix(X)) {
    stop("X is not a matrix. Data should be in matrix format.")
  }
  if (length(batch_vec) != nrow(X)) {
    stop("The number of rows in X and the number of batch labels are not equal.")
  }

  N <- nrow(X)
  P <- ncol(X)

  if (is_mixed) {
    if (is.null(column_type)) {
      stop("column_type must be supplied when type = 'MVN_MIXED'.")
    }
    if (length(column_type) != P) {
      stop("column_type must have one entry per column of X (length P).")
    }
  }

  if (!any(batch_vec == 0)) {
    batch_vec <- as.numeric(as.factor(batch_vec)) - 1
  }
  B <- length(unique(batch_vec))

  if (is.null(concentration)) {
    concentration <- rep(alpha, K)
  }

  # Empirical-Bayes hyperparameters. Mirrors mvnSampler::mvnSampler() /
  # mvnSamplerSeparationStrategy::mvnSamplerSeparationStrategy() exactly
  # (kappa is a fixed constant in both, never user-configurable - if that
  # ever changes there, it must change here too), including their
  # column-mean imputation of any missing entry for this one-off
  # calculation only (see imputeColumnMeans() in
  # src/genericFunctions.cpp) - this function never fits a model to X, so
  # there is no per-sweep augmentation step to mirror beyond that.
  kappa <- 0.01
  nu <- P + 2
  X_imputed <- X
  if (anyNA(X_imputed)) {
    col_means <- colMeans(X_imputed, na.rm = TRUE)
    for (p in seq_len(P)) {
      na_idx <- is.na(X_imputed[, p])
      X_imputed[na_idx, p] <- col_means[p]
    }
  }
  xi <- colMeans(X_imputed)
  global_cov <- stats::cov(X_imputed)
  iw_scale <- global_cov / K^(2 / P)
  delta_2 <- mean(diag(global_cov))
  lambda_2 <- m_scale
  # m(p, b) ~ N(batch_shift_prior_mean, delta_2 * lambda_2), i.e. standard
  # deviation sqrt(delta_2 * lambda_2) - see sampleMPrior() in
  # src/mvnSampler.cpp/mvnSamplerSeparationStrategy.cpp:
  # batch_shift_prior_precision there is a PRECISION (1 / (delta_2 *
  # lambda_2)), so the standard deviation multiplying the standard normal
  # draw is its inverse square root, not the precision (or its reciprocal)
  # directly.
  batch_shift_prior_scale <- sqrt(delta_2 * lambda_2)
  S_loc <- 1.0

  draw_cov <- function() {
    if (type == "MVN_LKJ" || is_mixed) {
      # Fixed (data-independent) priors: R_k ~ LKJ(eta), log(sigma_{k,p}) ~
      # N(beta, xi) with beta/xi the same fixed constants as
      # mvnSamplerSeparationStrategy.h (beta = 0.5 * log(0.72), xi = 1.0).
      R_k <- sampleLKJCorrelationMatrix(P, eta)
      sigma_k <- exp(stats::rnorm(P, mean = 0.5 * log(0.72), sd = 1.0))
      # Binary columns' marginal SD is fixed at 1 (the standard probit
      # identification device - see sampleCovPrior() in
      # src/mvnSamplerMixed.cpp); correlations are still free.
      if (is_mixed) {
        sigma_k[column_type == 1] <- 1.0
      }
      diag(sigma_k, P) %*% R_k %*% diag(sigma_k, P)
    } else {
      # Inverse-Wishart(iw_scale, nu).
      solve(stats::rWishart(1, nu, solve(iw_scale))[, , 1])
    }
  }

  simulate_one <- function() {
    # Prior cluster weights: w ~ Dirichlet(concentration), via the same
    # Gamma-normalisation construction sampler::updateWeights() itself uses
    # when every N_k is zero (i.e. before any data is assigned).
    g <- stats::rgamma(K, shape = concentration, rate = 1)
    w <- g / sum(g)

    mu <- matrix(0, P, K)
    cov <- array(0, dim = c(P, P, K))
    t_df <- rep(NA_real_, K)
    for (k in seq_len(K)) {
      cov[, , k] <- draw_cov()
      mu[, k] <- as.numeric(.mvtnormCholRnorm(xi, cov[, , k] / kappa))
      if (type == "MVT") {
        t_df[k] <- t_df_loc + stats::rgamma(1, shape = t_df_shape, rate = t_df_rate)
      }
    }

    m <- matrix(0, P, B)
    S <- matrix(0, P, B)
    for (b in seq_len(B)) {
      m[, b] <- stats::rnorm(P, mean = 0, sd = 1) * batch_shift_prior_scale
      S[, b] <- S_loc + 1 / stats::rgamma(P, shape = rho, rate = theta)
      # Binary columns' batch scale is fixed at 1 too (sampleSPrior() in
      # src/mvnSamplerMixed.cpp) - the same identification device as sigma.
      if (is_mixed) {
        S[column_type == 1, b] <- 1.0
      }
    }

    labels <- sample(seq_len(K), N, replace = TRUE, prob = w) - 1L

    X_sim <- matrix(NA_real_, N, P)
    for (n in seq_len(N)) {
      k <- labels[n] + 1L
      b <- batch_vec[n] + 1L
      mean_sum <- mu[, k] + m[, b]
      # matrix(..., nrow = P): cov[, , k] collapses to a bare scalar
      # (breaking diag<-, which needs an actual matrix) whenever P == 1 - a
      # perfectly ordinary univariate dataset.
      cov_comb <- matrix(cov[, , k], nrow = P)
      diag(cov_comb) <- diag(cov_comb) * S[, b]

      if (type == "MVT") {
        z <- .mvtnormCholRnorm(rep(0, P), cov_comb)
        w_chisq <- stats::rchisq(1, df = t_df[k]) / t_df[k]
        z <- mean_sum + as.numeric(z) / sqrt(w_chisq)
      } else {
        z <- mean_sum + as.numeric(.mvtnormCholRnorm(rep(0, P), cov_comb))
      }

      # Binary/probit columns: threshold the latent draw at 0 (Albert &
      # Chib, 1993) to get the replicated {0, 1} observation; continuous
      # columns keep the latent draw as-is.
      if (is_mixed) {
        z[column_type == 1] <- as.numeric(z[column_type == 1] > 0)
      }
      X_sim[n, ] <- z
    }

    list(
      X = X_sim,
      labels = labels,
      params = list(mu = mu, cov = cov, m = m, S = S, w = w, t_df = t_df)
    )
  }

  replicate(n_datasets, simulate_one(), simplify = FALSE)
}

#' @title Simulate from the posterior predictive distribution
#' @description For a sample of retained MCMC iterations, simulates a
#' replicate dataset from that iteration's drawn parameters and cluster
#' allocations - a posterior predictive check: "does the fitted model
#' generate datasets that look like the data we actually observed?" (Gelman
#' et al., 2013, \emph{BDA3}, Section 6.3; McElreath, 2020, \emph{Statistical
#' Rethinking}, 2nd ed., Section 3.3.2 for the single-parameter intuition,
#' Ch. 12/14 for hierarchical/multilevel posterior predictive checks in the
#' same spirit as this package's batch structure; Gelman et al., 2020,
#' "Bayesian Workflow", Section 2.2). Unlike \code{\link{simulatePriorPredictive}},
#' this uses parameter values already stored in \code{mcmc_output} directly
#' - no prior/empirical-Bayes formula is re-derived here, so there is no risk
#' of drifting out of sync with the sampler's own math.
#'
#' Every model type returned by \code{batchSemiSupervisedMixtureModel()} is
#' supported, via the stored \code{mean_sum}/\code{cov_comb} (and
#' \code{t_df} for MVT) arrays: \code{"MVN"}, \code{"MVT"}, \code{"MVN_LKJ"}
#' directly, and \code{"MVN_MIXED"} by additionally thresholding
#' binary/probit-linked columns (\code{column_type == 1}) at 0 (Albert &
#' Chib, 1993) after the same latent Gaussian draw - see
#' \code{\link{simulatePriorPredictive}}'s documentation for the same device
#' and why missing/censored cells are not specially handled here (every
#' replicate is fully observed; \code{\link{plotPredictiveCheck}}'s
#' \code{censor_code} argument handles the comparison side instead).
#' @param mcmc_output Output of \code{\link{batchSemiSupervisedMixtureModel}}/
#' \code{\link{runBatchMix}}.
#' @param batch_vec The batch label used to fit \code{mcmc_output} (0- or
#' 1-indexed; matched against \code{mcmc_output$B} either way).
#' @param column_type Only used if \code{mcmc_output$type == "MVN_MIXED"}: a
#' P-vector, \code{0} for a continuous column, \code{1} for a binary/probit
#' column - the same vector originally passed to
#' \code{\link{batchSemiSupervisedMixtureModel}}. Required for that type
#' (not stored in \code{mcmc_output} itself).
#' @param n_draws Number of posterior iterations to replicate from (sampled
#' without replacement from the retained, post-burn-in iterations).
#' @param burn Number of *original* MCMC iterations (not thinned samples) to
#' discard as burn-in before drawing from the remaining ones.
#' @param seed Optional seed for reproducible draw selection.
#' @return A list of length \code{n_draws}; each element is a list with
#' \code{X} (the simulated, fully-observed N x P replicate) and
#' \code{iteration} (the index, into the thinned/retained samples, the
#' replicate was drawn from).
#' @seealso \code{\link{simulatePriorPredictive}},
#' \code{\link{plotPredictiveCheck}}
#' @export
simulatePosteriorPredictive <- function(mcmc_output,
                                        batch_vec,
                                        column_type = NULL,
                                        n_draws = 50,
                                        burn = 0,
                                        seed = NULL) {
  is_mixed <- isTRUE(mcmc_output$type == "MVN_MIXED")
  P <- mcmc_output$P
  if (is_mixed) {
    if (is.null(column_type)) {
      stop("column_type must be supplied when mcmc_output$type == 'MVN_MIXED'.")
    }
    if (length(column_type) != P) {
      stop("column_type must have one entry per column (length mcmc_output$P).")
    }
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_iter <- mcmc_output$n_iter
  thin <- mcmc_output$thin
  n_saved <- nrow(mcmc_output$samples)
  K <- mcmc_output$K_max
  B <- mcmc_output$B
  is_mvt <- isTRUE(mcmc_output$type == "MVT")

  if (!any(batch_vec == 0)) {
    batch_vec <- as.numeric(as.factor(batch_vec)) - 1
  }
  N <- length(batch_vec)

  first_retained <- max(1L, ceiling((burn + thin) / thin))
  eligible <- seq(first_retained, n_saved)
  if (length(eligible) == 0) {
    stop("No retained iterations remain after applying `burn`.")
  }
  draw_idx <- sample(eligible, size = min(n_draws, length(eligible)), replace = FALSE)

  simulate_one <- function(r) {
    labels_r <- mcmc_output$samples[r, ]
    # matrix(..., nrow = P), not a bare [, , r] drop: indexing a single
    # slice off a 3D array silently collapses to a 1D vector whenever the
    # remaining (K * B) dimension is 1 too (e.g. K_max = B = 1), breaking
    # the column indexing below.
    mean_sum_r <- matrix(mcmc_output$mean_sum[, , r], nrow = P) # P x (K * B)
    cov_comb_r <- matrix(mcmc_output$cov_comb[, , r], nrow = P) # P x (P * K * B)

    X_sim <- matrix(NA_real_, N, P)
    for (n in seq_len(N)) {
      k <- labels_r[n]
      b <- batch_vec[n]
      kb <- k * B + b

      mean_sum_kb <- mean_sum_r[, kb + 1]
      cov_comb_kb <- cov_comb_r[, (kb * P + 1):(kb * P + P)]

      if (is_mvt) {
        t_df_k <- mcmc_output$t_df[r, k + 1]
        z <- .mvtnormCholRnorm(rep(0, P), cov_comb_kb)
        w_chisq <- stats::rchisq(1, df = t_df_k) / t_df_k
        z <- mean_sum_kb + as.numeric(z) / sqrt(w_chisq)
      } else {
        z <- mean_sum_kb + as.numeric(.mvtnormCholRnorm(rep(0, P), cov_comb_kb))
      }

      # Binary/probit columns: threshold the latent draw at 0 (Albert &
      # Chib, 1993); continuous columns keep the latent draw as-is.
      if (is_mixed) {
        z[column_type == 1] <- as.numeric(z[column_type == 1] > 0)
      }
      X_sim[n, ] <- z
    }

    list(X = X_sim, iteration = r)
  }

  lapply(draw_idx, simulate_one)
}

# Draws one row from N(mean, cov) via a Cholesky factor, tolerating a
# covariance matrix that is only positive semi-definite (a Metropolis
# proposal or an early-chain draw can occasionally leave cov_comb on the
# boundary) by falling back to Higham's nearest-PD-matrix correction rather
# than erroring. Kept dependency-free (no mvtnorm/MASS import) since this is
# the only place such a draw is needed.
.mvtnormCholRnorm <- function(mean, cov) {
  P <- length(mean)
  ch <- tryCatch(
    chol(cov),
    error = function(e) {
      eig <- eigen((cov + t(cov)) / 2, symmetric = TRUE)
      vals <- pmax(eig$values, 1e-8)
      chol(eig$vectors %*% diag(vals, P) %*% t(eig$vectors))
    }
  )
  mean + as.numeric(stats::rnorm(P) %*% ch)
}
