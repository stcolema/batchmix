#!/usr/bin/Rscript
#' @title Fit a batch mixture model
#' @description The package's main entry point: fits \code{n_chains}
#' independent chains of the batch mixture model and reports whether they
#' agree before handing anything back. Every Metropolis-Hastings proposal
#' window is auto-tuned by default (\code{auto_tune = TRUE}, via
#' Robbins-Monro diminishing adaptation over the first \code{n_burn}
#' iterations, frozen thereafter) - manually tuning \code{mu_proposal_window}
#' and friends is not required for typical use. Fitting more than one chain
#' (rather than calling \code{\link{runBatchMix}} once) is not optional
#' extra rigour - it is how non-convergence is detected at all: a single
#' chain can look perfectly stable on its own and still have converged to a
#' poor local optimum (mixture models are prone to this - see
#' \code{\link{generateInitialLabels}}'s k-means-seeded default, which
#' reduces but does not eliminate the risk). \code{n_chains >= 2}
#' automatically reports the rank-normalized, folded, split-Rhat and
#' bulk/tail ESS (Vehtari, Gelman, Simpson, Carpenter & Burkner, 2021,
#' \emph{Bayesian Analysis} 16(2)) on the permutation-invariant
#' complete-data log-likelihood trace, and identifies the chain with the
#' best (highest post-burn-in mean) BIC (see \code{\link{assessConvergence}}
#' and \code{\link{getBestChain}}).
#' @param n_chains Integer. Number of MCMC chains to run.
#' @param convergence_burn Number of iterations treated as burn-in when
#' assessing convergence across the chains just run (see
#' ``assessConvergence``); does not affect the returned chains themselves,
#' only the diagnostics attached to them. Defaults to half of ``n_iter``. Only
#' used if ``n_chains >= 2``.
#' @inheritParams runBatchMix
#' @param ... Further arguments passed to \code{\link{runBatchMix}} (e.g.
#' \code{include_interaction}, \code{batch_weight_prior} and their
#' associated options, or \code{r_proposal_window}/\code{sigma_proposal_window}/
#' \code{eta}/\code{column_type}/\code{censor_code} for \code{type}
#' \code{'MVN_LKJ'}/\code{'MVN_MIXED'} - see \code{?runBatchMix} for the full
#' list). Also accepts the deprecated \code{R} argument (renamed to
#' \code{n_iter}; still works, with a warning, for this release, when
#' passed positionally in its original slot or when every other argument is
#' also named - naming \code{R} while leaving later arguments positional is
#' not supported, since there is no way to bind two names to the same
#' argument slot).
#' @returns A list of named lists (one per chain, each the output of
#' ``runBatchMix``). If ``n_chains >= 2``, two attributes are attached to
#' the returned list: ``attr(., "convergence")`` (the output of
#' ``assessConvergence`` - rank-normalized split-Rhat, bulk/tail ESS, and
#' each chain's mean post-burn-in BIC) and ``attr(., "best_chain")`` (the
#' index of the chain with the highest mean post-burn-in BIC - see
#' ``getBestChain``). Passing the whole return value to
#' ``predictFromMultipleChains`` pools every chain regardless of this
#' attribute; use ``getBestChain()`` when a single representative chain
#' (rather than a pooled estimate) is wanted, e.g. to inspect or continue.
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
#' # Sampling parameters
#' n_iter <- 1000
#' thin <- 50
#' n_chains <- 4
#'
#' # MCMC samples
#' samples <- fitBatchMix(X, n_chains, n_iter, thin, batch_vec, "MVN",
#'   initial_labels = labels,
#'   fixed = fixed
#' )
#'
fitBatchMix <- function(X,
                        n_chains,
                        n_iter,
                        thin,
                        batch_vec,
                        type,
                        # -- problem specification --
                        K_max = NULL,
                        initial_labels = NULL,
                        fixed = NULL,
                        alpha = 1,
                        # -- MCMC control --
                        auto_tune = TRUE,
                        n_burn = NULL,
                        convergence_burn = NULL,
                        # -- proposal windows (only matter if auto_tune = FALSE) --
                        mu_proposal_window = 0.5**2,
                        cov_proposal_window = 0.002,
                        m_proposal_window = 0.3**2,
                        S_proposal_window = 0.01,
                        t_df_proposal_window = 0.015,
                        # -- prior hyperparameters --
                        m_scale = NULL,
                        rho = 3.0,
                        theta = 1.0,
                        # -- initial values (warm starts; default to prior draws) --
                        initial_class_means = NULL,
                        initial_class_covariance = NULL,
                        initial_batch_shift = NULL,
                        initial_batch_scale = NULL,
                        initial_class_df = NULL,
                        verbose = TRUE,
                        ...) {
  dots <- list(...)
  n_iter <- .resolveDeprecatedNIter(
    missing(n_iter), if (missing(n_iter)) NULL else n_iter, dots, "fitBatchMix"
  )
  # `R`, if present, has just been resolved into `n_iter` above - drop it so
  # it is not forwarded (and re-warned about) inside runBatchMix() too.
  dots$R <- NULL

  mcmc_lst <- vector("list", n_chains)

  mcmc_lst <- lapply(mcmc_lst, function(x) {
    do.call(runBatchMix, c(list(
      X,
      n_iter,
      thin,
      batch_vec,
      type,
      K_max = K_max,
      initial_labels = initial_labels,
      fixed = fixed,
      alpha = alpha,
      auto_tune = auto_tune,
      n_burn = n_burn,
      verbose = verbose,
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
      initial_class_df = initial_class_df
    ), dots))
  })

  # Record chain number
  for (ii in seq(n_chains)) {
    mcmc_lst[[ii]]$Chain <- ii
  }

  if (n_chains >= 2) {
    convergence <- tryCatch(
      assessConvergence(mcmc_lst, burn = convergence_burn),
      error = function(e) {
        if (verbose) {
          warning(paste0("Could not assess convergence across chains: ", conditionMessage(e)))
        }
        NULL
      }
    )

    if (!is.null(convergence)) {
      attr(mcmc_lst, "convergence") <- convergence
      attr(mcmc_lst, "best_chain") <- convergence$best_chain

      if (verbose) {
        message(sprintf(
          paste0(
            "\n%d chains run. Convergence (rank-normalized split-Rhat on %s, ",
            "post burn-in): Rhat = %.3f (bulk = %.3f, tail = %.3f), ",
            "ESS bulk = %.0f, ESS tail = %.0f.%s\nBest chain by mean BIC: #%d ",
            "(use getBestChain() to extract it).\n"
          ),
          n_chains, convergence$statistic,
          convergence$rhat, convergence$rhat_bulk, convergence$rhat_tail,
          convergence$ess_bulk, convergence$ess_tail,
          if (convergence$rhat > 1.01) {
            paste0(
              " Rhat > 1.01: these chains have likely NOT converged/mixed. A longer ",
              "chain sometimes helps, but persistent disagreement between chains ",
              "(especially unsupervised, well-separated clusters) is often a sign ",
              "that one or more chains are stuck in a different, worse-fitting local ",
              "optimum rather than merely under-run - re-fitting with more chains ",
              "(more chances at least one finds the better mode) is usually more ",
              "effective than simply lengthening the run."
            )
          } else {
            ""
          },
          convergence$best_chain
        ))
      }
    }
  }

  mcmc_lst
}
