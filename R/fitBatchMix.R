#!/usr/bin/Rscript
#' @title Fit a batch mixture model
#' @description The package's main entry point: fits \code{n_chains}
#' independent chains of the batch mixture model and reports whether they
#' agree before handing anything back. Every Metropolis-Hastings proposal
#' window (bundled into \code{control}, see \code{\link{batchmixControl}})
#' is auto-tuned by default (\code{control$auto_tune = TRUE}, via
#' Robbins-Monro diminishing adaptation over the first \code{control$n_burn}
#' iterations, frozen thereafter) - manually tuning proposal windows is not
#' required for typical use. Fitting more than one chain
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
#' @param auto_tune,n_burn,mu_proposal_window,cov_proposal_window,m_proposal_window,S_proposal_window,t_df_proposal_window
#' \strong{Deprecated}: pass these inside \code{control = batchmixControl(...)}
#' instead (e.g. \code{control = batchmixControl(auto_tune = FALSE)} rather
#' than \code{auto_tune = FALSE}). Still work this release (with a
#' warning); see \code{?batchmixControl} for what each one does. If both
#' \code{control} and one of these are supplied, \code{control} wins and
#' the deprecated argument is ignored (with a warning). (Repeated here,
#' rather than left to \code{@@inheritParams runBatchMix} above, because
#' these are also \code{fitBatchMix()}'s own formal arguments, not only
#' \code{runBatchMix()}'s, and roxygen2's \code{@@inheritParams} does not
#' resolve a comma-grouped \code{@@param} tag from another function.)
#' @param ... Further arguments passed to \code{\link{runBatchMix}} (e.g.
#' \code{include_interaction}, \code{batch_weight_prior} and their
#' associated options, or \code{eta}/\code{column_type}/\code{censor_code}
#' for \code{type} \code{'MVN_LKJ'}/\code{'MVN_MIXED'} - see
#' \code{?runBatchMix} for the full list). The MVN_LKJ/MVN_MIXED-only
#' proposal windows (\code{r_proposal_window}, \code{sigma_proposal_window},
#' \code{gamma_proposal_window}, \code{eta_proposal_window},
#' \code{gp_hyperparameter_proposal_window}) are also accepted here and
#' merged into \code{control} - passing them inside \code{control =
#' batchmixControl(...)} directly is preferred. Also accepts the deprecated
#' \code{R} argument (renamed to \code{n_iter}; still works, with a
#' warning, for this release, when passed positionally in its original slot
#' or when every other argument is also named - naming \code{R} while
#' leaving later arguments positional is not supported, since there is no
#' way to bind two names to the same argument slot).
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
                        control = batchmixControl(),
                        auto_tune = TRUE,
                        n_burn = NULL,
                        convergence_burn = NULL,
                        # -- proposal windows (deprecated - use `control` instead; only
                        # matter if auto_tune = FALSE) --
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
  # strict_dots = FALSE: unlike runBatchMix()/batchSemiSupervisedMixtureModel()
  # (whose `...` accepts only the deprecated `R`), fitBatchMix()'s `...` is
  # documented to forward arbitrary extra arguments on to runBatchMix() (e.g.
  # `include_interaction`, `batch_weight_prior`, `eta`, `column_type`) - only
  # `R` should be resolved/consumed here, everything else stays in `dots` to
  # be forwarded below.
  n_iter <- .resolveDeprecatedNIter(
    missing(n_iter), if (missing(n_iter)) NULL else n_iter, dots, "fitBatchMix",
    strict_dots = FALSE
  )
  # `R`, if present, has just been resolved into `n_iter` above - drop it so
  # it is not forwarded (and re-warned about) inside runBatchMix() too.
  dots$R <- NULL

  # `control` bundles every sampler-tuning argument (see ?batchmixControl).
  # mu_/cov_/m_/S_/t_df_proposal_window/auto_tune/n_burn are deprecated
  # formals here; r_/sigma_/gamma_/eta_/gp_hyperparameter_proposal_window
  # (MVN_LKJ/MVN_MIXED-only) are not formals of fitBatchMix() and instead
  # arrive via `...` - pulled out of `dots` here so they merge into
  # `control` too, rather than being forwarded raw (which would otherwise
  # silently override the resolved `control` once it reaches runBatchMix()).
  deprecated_control_args <- list()
  if (!missing(auto_tune)) deprecated_control_args$auto_tune <- auto_tune
  if (!missing(n_burn)) deprecated_control_args$n_burn <- n_burn
  if (!missing(mu_proposal_window)) deprecated_control_args$mu_proposal_window <- mu_proposal_window
  if (!missing(cov_proposal_window)) deprecated_control_args$cov_proposal_window <- cov_proposal_window
  if (!missing(m_proposal_window)) deprecated_control_args$m_proposal_window <- m_proposal_window
  if (!missing(S_proposal_window)) deprecated_control_args$S_proposal_window <- S_proposal_window
  if (!missing(t_df_proposal_window)) deprecated_control_args$t_df_proposal_window <- t_df_proposal_window

  dots_control_fields <- c(
    "r_proposal_window", "sigma_proposal_window", "gamma_proposal_window",
    "eta_proposal_window", "gp_hyperparameter_proposal_window"
  )
  for (field in intersect(dots_control_fields, names(dots))) {
    deprecated_control_args[[field]] <- dots[[field]]
    dots[[field]] <- NULL
  }

  control <- .resolveControlArgs(missing(control), control, deprecated_control_args, "fitBatchMix")

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
      control = control,
      verbose = verbose,
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
        # format.batchmix_convergence() (R/batchmixFitMethods.R) owns this
        # wording - also used by print.batchmix_convergence()/
        # print.batchmix_fit_list(), so the two never drift apart.
        message("\n", format(convergence, n_chains = n_chains))
      }
    }
  }

  # A thin S3 wrapper (still a plain list of batchmix_fit chains - every
  # existing `$`/`[[` access keeps working unchanged) so the console shows
  # a short report instead of dumping every chain's sampled arrays - see
  # R/batchmixFitMethods.R.
  class(mcmc_lst) <- c("batchmix_fit_list", class(mcmc_lst))

  mcmc_lst
}
