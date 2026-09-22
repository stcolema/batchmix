#!/usr/bin/Rscript
# S3 print()/summary()/format() methods for the classes fitBatchMix()/
# runBatchMix()/batchSemiSupervisedMixtureModel() (and the functions built
# on them - continueChain(s)(), processMCMCChain(s)()) attach to their
# return value: "batchmix_fit" (one chain), "batchmix_fit_list" (several
# chains) and "batchmix_convergence" (assessConvergence()'s return value).
#
# Every one of these objects is still a plain list underneath (a chain is
# still accessed with mcmc_out$mu, a chain list with mcmc_lst[[i]]) - the
# class only replaces the console's default behaviour of dumping every
# sampled array when a fit is printed, mclust/stanfit-style. See
# vignette("batchmix_workflow") for the full return-value shape.

# Mean acceptance rate per parameter family (Mu, Sigma/R/sigma, m, S, nu,
# gamma, eta, gp_hyperparameter - see collectAcceptanceRates()), collapsing
# the per-cluster/per-batch columns collectAcceptanceRates() returns (e.g.
# "Mu_1".."Mu_K") down to one number per family for a compact print().
.acceptanceRateFamilyMeans <- function(rates_df) {
  nm <- names(rates_df)
  if (length(nm) == 0) {
    return(NULL)
  }
  family <- sub("_[0-9]+$", "", nm)
  fam_levels <- unique(family)
  means <- vapply(fam_levels, function(f) {
    mean(as.matrix(rates_df[, family == f, drop = FALSE]), na.rm = TRUE)
  }, numeric(1))
  names(means) <- fam_levels
  means
}

#' @title Print a batchmix fit
#' @description Prints a short report (density type, dimensions, MCMC
#' settings, mean acceptance rates and, if attached, convergence
#' diagnostics) instead of dumping the fit's raw sampled arrays -
#' \code{summary()} gives more detail; \code{x} remains a plain list
#' underneath (\code{x$mu} etc. keep working) - the class only changes how
#' it prints.
#' @param x Output of \code{\link{runBatchMix}},
#' \code{\link{batchSemiSupervisedMixtureModel}}, \code{\link{continueChain}}
#' or \code{\link{processMCMCChain}} (one chain).
#' @param ... Unused; present for S3 method consistency.
#' @return \code{x}, invisibly.
#' @export
print.batchmix_fit <- function(x, ...) {
  supervision <- if (isTRUE(x$Semisupervised)) "semi-supervised" else "unsupervised"
  n_retained <- tryCatch(floor(x$n_iter / x$thin), error = function(e) NA_integer_)

  cat(sprintf("<batchmix fit> type = %s, %s\n", x$type, supervision))
  cat(sprintf("  N = %d, P = %d, K_max = %d, B = %d\n", x$N, x$P, x$K_max, x$B))
  cat(sprintf(
    "  n_iter = %d, thin = %d (%d samples retained), auto_tune = %s\n",
    x$n_iter, x$thin, n_retained, isTRUE(x$auto_tune)
  ))

  rates <- tryCatch(collectAcceptanceRates(x), error = function(e) NULL)
  fam_means <- if (!is.null(rates)) .acceptanceRateFamilyMeans(rates) else NULL
  if (!is.null(fam_means)) {
    cat(
      "  Mean acceptance rates:",
      paste(sprintf("%s = %.2f", names(fam_means), fam_means), collapse = ", "),
      "\n"
    )
  }

  convergence <- attr(x, "convergence")
  if (!is.null(convergence)) {
    cat("\n")
    print(convergence)
  }

  cat("\nUse summary() for more detail.\n")
  invisible(x)
}

#' @title Summarise a batchmix fit
#' @description Collects dimensions, MCMC settings, BIC/likelihood and mean
#' acceptance rates for a single fitted chain into a compact object; print
#' the result for a formatted report.
#' @param object Output of \code{\link{runBatchMix}},
#' \code{\link{batchSemiSupervisedMixtureModel}}, \code{\link{continueChain}}
#' or \code{\link{processMCMCChain}} (one chain).
#' @param ... Unused; present for S3 method consistency.
#' @return An object of class \code{summary.batchmix_fit}.
#' @export
summary.batchmix_fit <- function(object, ...) {
  rates <- tryCatch(collectAcceptanceRates(object), error = function(e) NULL)
  fam_means <- if (!is.null(rates)) .acceptanceRateFamilyMeans(rates) else NULL

  bic <- object$BIC
  cll <- object$complete_likelihood
  oll <- object$observed_likelihood

  structure(
    list(
      type = object$type,
      Semisupervised = isTRUE(object$Semisupervised),
      N = object$N, P = object$P, K_max = object$K_max, B = object$B,
      n_iter = object$n_iter, thin = object$thin, n_burn = object$n_burn,
      auto_tune = isTRUE(object$auto_tune),
      bic_mean = if (!is.null(bic)) mean(bic) else NA_real_,
      bic_range = if (!is.null(bic)) range(bic) else c(NA_real_, NA_real_),
      complete_likelihood_mean = if (!is.null(cll)) mean(cll) else NA_real_,
      observed_likelihood_mean = if (!is.null(oll)) mean(oll) else NA_real_,
      acceptance_rates = fam_means,
      convergence = attr(object, "convergence")
    ),
    class = "summary.batchmix_fit"
  )
}

#' @rdname summary.batchmix_fit
#' @param x Output of \code{summary.batchmix_fit}.
#' @export
print.summary.batchmix_fit <- function(x, ...) {
  cat(sprintf(
    "<batchmix fit summary> type = %s, %s\n",
    x$type, if (x$Semisupervised) "semi-supervised" else "unsupervised"
  ))
  cat(sprintf("  N = %d, P = %d, K_max = %d, B = %d\n", x$N, x$P, x$K_max, x$B))
  cat(sprintf(
    "  n_iter = %d, thin = %d, n_burn = %s, auto_tune = %s\n",
    x$n_iter, x$thin, if (is.null(x$n_burn)) "NA" else x$n_burn, x$auto_tune
  ))
  cat(sprintf(
    "  BIC: mean = %.2f, range = [%.2f, %.2f]\n",
    x$bic_mean, x$bic_range[1], x$bic_range[2]
  ))
  cat(sprintf(
    "  Complete-data log-likelihood (mean) = %.2f, observed log-likelihood (mean) = %.2f\n",
    x$complete_likelihood_mean, x$observed_likelihood_mean
  ))
  if (!is.null(x$acceptance_rates)) {
    cat("  Mean acceptance rates by parameter family:\n")
    for (nm in names(x$acceptance_rates)) {
      cat(sprintf("    %-18s %.2f\n", nm, x$acceptance_rates[[nm]]))
    }
  }
  if (!is.null(x$convergence)) {
    cat("\n")
    print(x$convergence)
  }
  invisible(x)
}

#' @title Print a list of batchmix fits
#' @description Prints a short report (chain count, dimensions, MCMC
#' settings and, if computed, convergence diagnostics) instead of dumping
#' every chain's raw sampled arrays - \code{summary()} gives per-chain
#' detail; \code{x} remains a plain list of chains underneath
#' (\code{x[[i]]$mu} etc. keep working) - the class only changes how it
#' prints.
#' @param x Output of \code{\link{fitBatchMix}}, \code{\link{continueChains}}
#' or \code{\link{processMCMCChains}}.
#' @param ... Unused; present for S3 method consistency.
#' @return \code{x}, invisibly.
#' @export
print.batchmix_fit_list <- function(x, ...) {
  n_chains <- length(x)
  first <- x[[1]]

  cat(sprintf(
    "<batchmix fit list> %d chain%s, type = %s, %s\n",
    n_chains, if (n_chains == 1) "" else "s", first$type,
    if (isTRUE(first$Semisupervised)) "semi-supervised" else "unsupervised"
  ))
  cat(sprintf("  N = %d, P = %d, K_max = %d, B = %d\n", first$N, first$P, first$K_max, first$B))
  cat(sprintf("  n_iter = %d, thin = %d\n", first$n_iter, first$thin))

  convergence <- attr(x, "convergence")
  if (!is.null(convergence)) {
    cat("\n")
    print(convergence, n_chains = n_chains)
  } else if (n_chains >= 2) {
    cat("\n  Convergence not computed (assessConvergence() failed, or was not run).\n")
  }

  cat("\nUse summary() for more detail, or getBestChain()/processMCMCChains() to continue.\n")
  invisible(x)
}

#' @title Summarise a list of batchmix fits
#' @description Collects a \code{\link{summary.batchmix_fit}} for each
#' chain, plus the attached convergence diagnostics and best-chain index (if
#' any); print the result for a formatted report.
#' @param object Output of \code{\link{fitBatchMix}},
#' \code{\link{continueChains}} or \code{\link{processMCMCChains}}.
#' @param ... Unused; present for S3 method consistency.
#' @return An object of class \code{summary.batchmix_fit_list}.
#' @export
summary.batchmix_fit_list <- function(object, ...) {
  structure(
    list(
      n_chains = length(object),
      per_chain = lapply(object, summary),
      convergence = attr(object, "convergence"),
      best_chain = attr(object, "best_chain")
    ),
    class = "summary.batchmix_fit_list"
  )
}

#' @rdname summary.batchmix_fit_list
#' @param x Output of \code{summary.batchmix_fit_list}.
#' @export
print.summary.batchmix_fit_list <- function(x, ...) {
  cat(sprintf("<batchmix fit list summary> %d chains\n", x$n_chains))

  bic_means <- vapply(x$per_chain, function(s) s$bic_mean, numeric(1))
  cat(
    "  Per-chain mean BIC:",
    paste(sprintf("#%d = %.2f", seq_along(bic_means), bic_means), collapse = ", "),
    "\n"
  )
  if (!is.null(x$best_chain)) {
    cat(sprintf("  Best chain by mean BIC: #%d\n", x$best_chain))
  }

  if (!is.null(x$convergence)) {
    cat("\n")
    print(x$convergence, n_chains = x$n_chains)
  }
  invisible(x)
}

#' @title Format an assessConvergence() result
#' @description Builds the Rhat/ESS/best-chain report text shared by
#' \code{print.batchmix_convergence()}, \code{print.batchmix_fit()},
#' \code{print.batchmix_fit_list()} and \code{\link{fitBatchMix}}'s
#' automatic post-fit message, so the wording cannot drift out of sync
#' between them.
#' @param x Output of \code{\link{assessConvergence}}.
#' @param n_chains If given, prefixed as "\code{n_chains} chains run." -
#' used right after fitting, omitted when reporting on an
#' already-fitted/attached object.
#' @param ... Unused; present for S3 method consistency.
#' @return A length-one character string.
#' @export
format.batchmix_convergence <- function(x, n_chains = NULL, ...) {
  chain_note <- if (!is.null(n_chains)) sprintf("%d chains run. ", n_chains) else ""
  warn <- if (x$rhat > 1.01) {
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
  }

  sprintf(
    paste0(
      "%sConvergence (rank-normalized split-Rhat on %s, post burn-in): ",
      "Rhat = %.3f (bulk = %.3f, tail = %.3f), ESS bulk = %.0f, ESS tail = %.0f.%s\n",
      "Best chain by mean BIC: #%d (use getBestChain() to extract it)."
    ),
    chain_note, x$statistic,
    x$rhat, x$rhat_bulk, x$rhat_tail, x$ess_bulk, x$ess_tail, warn,
    x$best_chain
  )
}

#' @title Print an assessConvergence() result
#' @description Prints the Rhat/ESS/best-chain report built by
#' \code{\link{format.batchmix_convergence}}.
#' @param x Output of \code{\link{assessConvergence}}.
#' @param ... Passed on to \code{\link{format.batchmix_convergence}} (e.g.
#' \code{n_chains}).
#' @return \code{x}, invisibly.
#' @export
print.batchmix_convergence <- function(x, ...) {
  cat(format(x, ...), "\n")
  invisible(x)
}
