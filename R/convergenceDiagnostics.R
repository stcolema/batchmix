#' @title Rank-normalized split-Rhat and effective sample size
#' @description Computes the rank-normalized, folded, split-\eqn{\hat{R}}
#' convergence diagnostic and the associated bulk/tail effective sample
#' sizes (ESS) for a set of parallel MCMC chains, following Vehtari, Gelman,
#' Simpson, Carpenter & Burkner (2021), "Rank-normalization, folding, and
#' localization: An improved \eqn{\hat{R}} for assessing convergence of
#' MCMC", \emph{Bayesian Analysis} 16(2), 667-718 - the diagnostic that
#' superseded the classical Gelman-Rubin \eqn{\hat{R}} (Gelman & Rubin,
#' 1992) precisely because the classical version can fail to detect
#' non-convergence in the tails of a distribution, or when chains have very
#' different (but individually stable) variances. This is what
#' \code{rstan}/\code{posterior}/\code{bayesplot} report by default as of
#' Stan 2.29+.
#'
#' Rank-normalizing (replacing each draw by its rank across all chains,
#' transformed to a standard normal quantile) before computing the
#' classical split-\eqn{\hat{R}} formula makes the diagnostic robust to
#' heavy tails and insensitive to the parameter's actual scale; "folding"
#' (repeating the same computation on \eqn{|\theta - \mathrm{median}(\theta)|}
#' and taking the worse of the two) additionally catches chains that agree
#' on location but disagree on scale. \eqn{\hat{R}} values comfortably above
#' 1.01 indicate chains that have not mixed (Vehtari et al., 2021,
#' recommend 1.01 as the threshold, tightened from the older folklore
#' value of 1.1).
#'
#' @param chains A numeric matrix, iterations (rows) x chains (columns), OR
#' a list of numeric vectors (one per chain, need not be equal length after
#' removing burn-in elsewhere - they are truncated to the shortest).
#' @return A list with \code{rhat} (the reported \eqn{\hat{R}}, i.e.
#' \code{max(rhat_bulk, rhat_tail)}), \code{rhat_bulk}, \code{rhat_tail},
#' \code{ess_bulk} and \code{ess_tail}.
#' @references Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. &
#' Burkner, P-C. (2021). Rank-normalization, folding, and localization: An
#' improved \eqn{\hat{R}} for assessing convergence of MCMC. \emph{Bayesian
#' Analysis}, 16(2), 667-718.
#' @export
#' @examples
#' set.seed(1)
#' # Three well-mixed chains targeting the same distribution
#' good_chains <- matrix(rnorm(3000), ncol = 3)
#' rankNormalizedRhat(good_chains)$rhat
#'
#' # Three chains that have not mixed (different means)
#' bad_chains <- cbind(rnorm(1000, 0), rnorm(1000, 2), rnorm(1000, -2))
#' rankNormalizedRhat(bad_chains)$rhat
rankNormalizedRhat <- function(chains) {
  if (is.list(chains)) {
    min_len <- min(vapply(chains, length, integer(1)))
    chains <- vapply(chains, function(x) utils::tail(x, min_len), numeric(min_len))
  }
  if (!is.matrix(chains)) {
    stop("`chains` must be a matrix (iterations x chains) or a list of numeric vectors.")
  }

  n <- nrow(chains)
  m <- ncol(chains)
  if (n < 4) {
    stop("Need at least 4 iterations per chain (after burn-in) to split.")
  }

  # Split each chain in half (split-Rhat): doubles the number of chains
  # compared, at half the length, so within-chain non-stationarity shows up
  # as a between-"chain" discrepancy too.
  half <- n %/% 2
  split_chains <- cbind(
    chains[seq_len(half), , drop = FALSE],
    chains[(n - half + 1):n, , drop = FALSE]
  )

  rhat_ess <- function(x) {
    n_i <- nrow(x)
    m_i <- ncol(x)

    z <- .rankNormalize(x)

    chain_means <- colMeans(z)
    grand_mean <- mean(chain_means)
    B <- n_i / (m_i - 1) * sum((chain_means - grand_mean)^2)
    W <- mean(apply(z, 2, stats::var))
    var_plus <- ((n_i - 1) / n_i) * W + B / n_i
    rhat <- sqrt(var_plus / W)

    ess <- .effectiveSampleSize(z, W, var_plus)

    list(rhat = rhat, ess = ess)
  }

  bulk <- rhat_ess(split_chains)
  folded <- abs(split_chains - stats::median(split_chains))
  tail_metric <- rhat_ess(folded)

  list(
    rhat = max(bulk$rhat, tail_metric$rhat),
    rhat_bulk = bulk$rhat,
    rhat_tail = tail_metric$rhat,
    ess_bulk = bulk$ess,
    ess_tail = tail_metric$ess
  )
}

# Average-rank transform to the standard normal quantile scale (Blom's
# formula, matching the 'posterior' package's default), applied across ALL
# chains jointly (rank-normalization must pool chains to be comparable).
.rankNormalize <- function(x) {
  r <- rank(x, ties.method = "average")
  z <- stats::qnorm((r - 3 / 8) / (length(r) - 1 / 4))
  matrix(z, nrow = nrow(x), ncol = ncol(x))
}

# Multi-chain effective sample size via the sum-of-autocorrelations
# estimator (Gelman et al., 2013, BDA3, Eq. 11.8), truncating the
# autocorrelation sum at Geyer's (1992) initial monotone sequence: sum
# consecutive PAIRS of lag-autocorrelations while their sum stays positive
# and non-increasing, discard the rest. This is the same construction
# 'rstan'/'posterior' use for bulk/tail ESS once the input has already been
# rank-normalized (bulk) or rank-normalized-after-folding (tail).
.effectiveSampleSize <- function(z, W, var_plus) {
  n_i <- nrow(z)
  m_i <- ncol(z)

  # Too few post-split iterations to estimate any autocorrelation structure
  # at all (needs at least a lag-2 pair) - fall back to the no-autocorrelation
  # assumption (tau_hat = 1) rather than erroring; assessConvergence()'s
  # Rhat is still meaningful here even when ESS can't be, and a short chain
  # should read as "low ESS", not crash.
  max_lag <- n_i - 1
  if (max_lag < 2) {
    return(m_i * n_i)
  }

  acov <- vapply(seq_len(m_i), function(j) {
    stats::acf(z[, j], lag.max = max_lag, plot = FALSE, type = "covariance")$acf[, 1, 1]
  }, numeric(n_i))

  rho_hat_t <- 1 - (W - rowMeans(acov)) / var_plus

  # Geyer's initial monotone sequence: pair up lags (1,2), (3,4), ...; sum of
  # each pair must be positive and the running sum non-increasing.
  lag_starts <- seq(2, max_lag, by = 2)
  paired_sums <- rho_hat_t[lag_starts] + rho_hat_t[pmin(lag_starts + 1, max_lag + 1)]
  paired_sums[is.na(paired_sums)] <- -Inf

  keep <- integer(0)
  running_min <- Inf
  for (i in seq_along(paired_sums)) {
    if (paired_sums[i] <= 0) break
    running_min <- min(running_min, paired_sums[i])
    keep <- c(keep, running_min)
  }

  tau_hat <- 1 + 2 * sum(keep)
  ess <- (m_i * n_i) / max(tau_hat, 1e-8)
  min(ess, m_i * n_i)
}

#' @title Assess MCMC convergence across chains and identify the best chain
#' @description The single entry point for the convergence-checking part of
#' the Bayesian workflow this package supports: computes the
#' rank-normalized split-\eqn{\hat{R}} and bulk/tail ESS (see
#' \code{\link{rankNormalizedRhat}}) for a set of chains from
#' \code{\link{runMCMCChains}}, and identifies which chain has the best
#' (highest, post-burn-in mean) BIC.
#'
#' \strong{Why the complete-data log-likelihood/BIC trace, not the raw
#' cluster parameters}: mixture models are only identified up to a
#' permutation of the cluster labels ("label switching") - cluster 1 in one
#' chain may correspond to cluster 2 in another, purely by chance, with no
#' bearing on convergence. Computing \eqn{\hat{R}} directly on the raw
#' \code{mu}/\code{cov}/\code{m}/\code{S} arrays would conflate genuine
#' non-convergence with harmless label permutation and can be badly
#' misleading. The complete-data log-likelihood (and BIC, a fixed
#' transformation of it) is invariant to label permutation, so it is the
#' safe default target for this diagnostic; see Stephens (2000), "Dealing
#' with label switching in mixture models", JRSS-B 62(4), and Celeux, Hurn
#' & Robert (2000), JASA 95(451), for the label-switching problem in
#' general, and \code{\link{minVI}} for point-estimate summarisation that
#' works around it via the posterior similarity matrix instead of the raw
#' labels.
#' @param mcmc_chains Output of \code{\link{runMCMCChains}}.
#' @param burn Number of iterations to discard as burn-in before computing
#' diagnostics. Defaults to half the chain length.
#' @param statistic Which permutation-invariant trace to assess
#' convergence on: \code{"complete_likelihood"} (default),
#' \code{"observed_likelihood"} or \code{"BIC"}.
#' @return A list with \code{rhat}, \code{ess_bulk}, \code{ess_tail} (from
#' \code{\link{rankNormalizedRhat}} on \code{statistic}), \code{chain_bic}
#' (a vector of each chain's post-burn-in mean BIC), \code{best_chain} (the
#' index of the chain with the highest mean BIC) and \code{statistic} (which
#' trace was used).
#' @export
assessConvergence <- function(mcmc_chains,
                              burn = NULL,
                              statistic = c("complete_likelihood", "observed_likelihood", "BIC")) {
  statistic <- match.arg(statistic)
  n_chains <- length(mcmc_chains)
  if (n_chains < 2) {
    stop("Need at least 2 chains to assess convergence.")
  }

  n_iter <- mcmc_chains[[1]]$n_iter
  thin <- mcmc_chains[[1]]$thin
  if (is.null(burn)) {
    burn <- floor(n_iter / 2)
  }
  first_retained <- max(1L, ceiling((burn + thin) / thin))

  traces <- lapply(mcmc_chains, function(ch) {
    n_saved <- length(ch[[statistic]])
    ch[[statistic]][min(first_retained, n_saved):n_saved]
  })

  diag <- rankNormalizedRhat(traces)

  chain_bic <- vapply(mcmc_chains, function(ch) {
    n_saved <- length(ch$BIC)
    mean(ch$BIC[min(first_retained, n_saved):n_saved])
  }, numeric(1))

  structure(
    list(
      rhat = diag$rhat,
      rhat_bulk = diag$rhat_bulk,
      rhat_tail = diag$rhat_tail,
      ess_bulk = diag$ess_bulk,
      ess_tail = diag$ess_tail,
      chain_bic = chain_bic,
      best_chain = which.max(chain_bic),
      statistic = statistic
    ),
    class = "batchmix_convergence"
  )
}

#' @title Extract the best chain from a set of MCMC chains
#' @description Returns the single chain identified as "best" (highest
#' post-burn-in mean BIC) by \code{\link{assessConvergence}} - a convenience
#' for any downstream step that needs one chain rather than a pooled set
#' (\code{\link{predictFromMultipleChains}} pools all chains instead, which
#' is usually preferable when every chain has converged; use this when you
#' specifically want a single representative chain, e.g. to inspect or
#' \code{\link{continueChain}}). If \code{mcmc_chains} was produced by
#' \code{\link{runMCMCChains}}, its attached \code{"best_chain"}/
#' \code{"convergence"} attributes are reused rather than recomputed;
#' otherwise \code{\link{assessConvergence}} is called on the fly.
#' @param mcmc_chains Output of \code{\link{runMCMCChains}}.
#' @param burn Only used if convergence has not already been attached/
#' computed; see \code{\link{assessConvergence}}.
#' @return The single best chain (a named list, as returned by
#' \code{\link{batchSemiSupervisedMixtureModel}}), with its convergence
#' diagnostics attached as \code{attr(., "convergence")}.
#' @export
getBestChain <- function(mcmc_chains, burn = NULL) {
  convergence <- attr(mcmc_chains, "convergence")
  best_chain <- attr(mcmc_chains, "best_chain")

  if (is.null(convergence) || is.null(best_chain)) {
    convergence <- assessConvergence(mcmc_chains, burn = burn)
    best_chain <- convergence$best_chain
  }

  out <- mcmc_chains[[best_chain]]
  attr(out, "convergence") <- convergence
  out
}
