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
#' \strong{Only 'MVN', 'MVT' and 'MVN_LKJ' are supported} - 'MVN_MIXED' is
#' not, since its prior predictive draw would additionally need to simulate
#' the probit/missingness/censoring layer, which is out of scope here.
#'
#' @param X Data matrix (items in rows), used only to derive the
#' empirical-Bayes prior location/scale (its own values are never reused
#' directly) - the same summary statistics
#' (\code{colMeans(X)}/\code{cov(X)}) that \code{batchSemiSupervisedMixtureModel()}
#' computes internally when \code{type} is 'MVN', 'MVT' or 'MVN_LKJ'.
#' @param batch_vec Observed batch label for each row of X.
#' @param K The number of clusters to simulate.
#' @param type One of 'MVN', 'MVT', 'MVN_LKJ' (see
#' \code{\link{batchSemiSupervisedMixtureModel}} for what each means).
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
#' 'MVN_LKJ'.
#' @param t_df_shape,t_df_rate,t_df_loc Hyperparameters of the (shifted
#' Gamma) prior on the cluster degrees of freedom; only used if \code{type}
#' is 'MVT'. Defaults match \code{mvtSampler}'s own defaults
#' (\code{psi = 2, chi = 0.1, t_loc = 2}).
#' @param n_datasets Number of independent prior predictive datasets to
#' draw.
#' @return A list of length \code{n_datasets}; each element is a list with
#' \code{X} (the simulated N x P data matrix), \code{labels} (the simulated
#' N-vector of prior cluster draws, 0-indexed) and \code{params} (the drawn
#' mu/cov/m/S/weights, and t_df if \code{type == "MVT"}).
#' @seealso \code{\link{simulatePosteriorPredictive}},
#' \code{\link{plotPredictiveCheck}}
#' @export
simulatePriorPredictive <- function(X,
                                    batch_vec,
                                    K,
                                    type = c("MVN", "MVT", "MVN_LKJ"),
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

  if (!is.matrix(X)) {
    stop("X is not a matrix. Data should be in matrix format.")
  }
  if (anyNA(X)) {
    stop(paste0(
      "simulatePriorPredictive() does not support missing data (X contains ",
      "NA/NaN) for type = '", type, "' - this mirrors the same limitation ",
      "'MVN_MIXED' alone lifts in batchSemiSupervisedMixtureModel()."
    ))
  }
  if (length(batch_vec) != nrow(X)) {
    stop("The number of rows in X and the number of batch labels are not equal.")
  }

  N <- nrow(X)
  P <- ncol(X)

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
  # ever changes there, it must change here too).
  kappa <- 0.01
  nu <- P + 2
  xi <- colMeans(X)
  global_cov <- stats::cov(X)
  iw_scale <- global_cov / K^(2 / P)
  delta_2 <- mean(diag(global_cov))
  lambda_2 <- m_scale
  # Matches the C++ literally: m(p, b) is drawn as
  # batch_shift_prior_mean + batch_shift_prior_precision * Z, Z ~ N(0, 1) -
  # i.e. despite the name, this quantity multiplies the standard normal
  # draw directly (see src/mvnSampler.cpp::sampleMPrior()).
  batch_shift_prior_scale <- 1 / (delta_2 * lambda_2)
  S_loc <- 1.0

  draw_cov <- function() {
    if (type == "MVN_LKJ") {
      # Fixed (data-independent) priors: R_k ~ LKJ(eta), log(sigma_{k,p}) ~
      # N(beta, xi) with beta/xi the same fixed constants as
      # mvnSamplerSeparationStrategy.h (beta = 0.5 * log(0.72), xi = 1.0).
      R_k <- sampleLKJCorrelationMatrix(P, eta)
      sigma_k <- exp(stats::rnorm(P, mean = 0.5 * log(0.72), sd = 1.0))
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
    }

    labels <- sample(seq_len(K), N, replace = TRUE, prob = w) - 1L

    X_sim <- matrix(NA_real_, N, P)
    for (n in seq_len(N)) {
      k <- labels[n] + 1L
      b <- batch_vec[n] + 1L
      mean_sum <- mu[, k] + m[, b]
      cov_comb <- cov[, , k]
      diag(cov_comb) <- diag(cov_comb) * S[, b]

      if (type == "MVT") {
        z <- .mvtnormCholRnorm(rep(0, P), cov_comb)
        w_chisq <- stats::rchisq(1, df = t_df[k]) / t_df[k]
        X_sim[n, ] <- mean_sum + as.numeric(z) / sqrt(w_chisq)
      } else {
        X_sim[n, ] <- mean_sum + as.numeric(.mvtnormCholRnorm(rep(0, P), cov_comb))
      }
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
#' supported for the continuous case (\code{"MVN"}, \code{"MVT"},
#' \code{"MVN_LKJ"}) via the stored \code{mean_sum}/\code{cov_comb} (and
#' \code{t_df} for MVT) arrays; \code{"MVN_MIXED"} is not currently
#' supported (its replicate would additionally need the probit/missingness/
#' censoring observation layer).
#' @param mcmc_output Output of \code{\link{batchSemiSupervisedMixtureModel}}/
#' \code{\link{runBatchMix}} (\code{type} one of 'MVN', 'MVT', 'MVN_LKJ').
#' @param batch_vec The batch label used to fit \code{mcmc_output} (0- or
#' 1-indexed; matched against \code{mcmc_output$B} either way).
#' @param n_draws Number of posterior iterations to replicate from (sampled
#' without replacement from the retained, post-burn-in iterations).
#' @param burn Number of *original* MCMC iterations (not thinned samples) to
#' discard as burn-in before drawing from the remaining ones.
#' @param seed Optional seed for reproducible draw selection.
#' @return A list of length \code{n_draws}; each element is a list with
#' \code{X} (the simulated N x P replicate) and \code{iteration} (the index,
#' into the thinned/retained samples, the replicate was drawn from).
#' @seealso \code{\link{simulatePriorPredictive}},
#' \code{\link{plotPredictiveCheck}}
#' @export
simulatePosteriorPredictive <- function(mcmc_output,
                                        batch_vec,
                                        n_draws = 50,
                                        burn = 0,
                                        seed = NULL) {
  if (isTRUE(mcmc_output$type == "MVN_MIXED")) {
    stop("simulatePosteriorPredictive() does not support type = 'MVN_MIXED' yet.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_iter <- mcmc_output$n_iter
  thin <- mcmc_output$thin
  n_saved <- nrow(mcmc_output$samples)
  P <- mcmc_output$P
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
    mean_sum_r <- mcmc_output$mean_sum[, , r] # P x (K * B)
    cov_comb_r <- mcmc_output$cov_comb[, , r] # P x (P * K * B)

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
        X_sim[n, ] <- mean_sum_kb + as.numeric(z) / sqrt(w_chisq)
      } else {
        X_sim[n, ] <- mean_sum_kb + as.numeric(.mvtnormCholRnorm(rep(0, P), cov_comb_kb))
      }
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
