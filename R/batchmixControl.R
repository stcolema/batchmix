#!/usr/bin/Rscript
#' @title Sampler-tuning control for a batchmix fit
#' @description Bundles the Metropolis-Hastings proposal windows and the
#' auto-tuning schedule - the arguments that control \emph{how hard the
#' sampler works to explore the posterior}, not what is being fitted - into
#' one object, the way \code{stats::glm(control = glm.control(...))} or
#' \code{lme4::glmer(control = glmerControl(...))} do. Every proposal window
#' is only a starting value for Robbins-Monro auto-tuning
#' (\code{auto_tune = TRUE}, the default - see \code{\link{runBatchMix}});
#' getting them exactly right is not required for typical use, which is
#' also why they live here rather than among \code{\link{runBatchMix}}'s
#' main arguments.
#'
#' Prior/hyperparameter choices that change \emph{what model is fitted}
#' (\code{rho}, \code{theta}, \code{m_scale}, \code{eta},
#' \code{a_gamma}/\code{b_gamma}, \code{batch_weight_prior} and its
#' associated GP/partial-pooling hyperparameters) are deliberately not part
#' of \code{control} - they stay as ordinary named arguments to
#' \code{\link{runBatchMix}}/\code{\link{fitBatchMix}}, so that changing
#' them is never mistaken for a purely-technical, sampler-internal tweak.
#' @param mu_proposal_window The proposal window for the cluster mean
#' proposal kernel. The proposal density is a Gaussian distribution, the
#' window is the variance. Making this smaller will normally increase the
#' acceptance rate.
#' @param cov_proposal_window The proposal window for the cluster
#' covariance proposal kernel when \code{type} is \code{'MVN'} or
#' \code{'MVT'}. The proposal density is a Wishart distribution, this
#' argument is the reciprocal of the degree of freedom.
#' @param r_proposal_window Only used if \code{type} is \code{'MVN_LKJ'} or
#' \code{'MVN_MIXED'}: the standard deviation of the (unconstrained-space)
#' Gaussian random walk proposal for the cluster correlation matrix R.
#' @param sigma_proposal_window Only used if \code{type} is
#' \code{'MVN_LKJ'} or \code{'MVN_MIXED'}: the proposal window for the
#' cluster marginal standard deviations.
#' @param m_proposal_window The proposal window for the batch mean proposal
#' kernel. The proposal density is a Gaussian distribution, the window is
#' the variance.
#' @param S_proposal_window The proposal window for the batch standard
#' deviation proposal kernel. The proposal density is a Gamma distribution,
#' this argument is the reciprocal of the rate.
#' @param t_df_proposal_window The proposal window for the degrees of
#' freedom for the multivariate t distribution (only used if \code{type} is
#' \code{'MVT'}).
#' @param gamma_proposal_window Proposal window (Gaussian random-walk SD)
#' for the interaction term; only used if \code{include_interaction} is
#' \code{TRUE}.
#' @param eta_proposal_window Proposal window for the batch-weight
#' Metropolis-Hastings update; only used if \code{batch_weight_prior} is
#' \code{"partial_pooling"} or \code{"gp"}.
#' @param gp_hyperparameter_proposal_window Proposal window for the GP
#' hyperparameter update; only used if \code{batch_weight_prior} is
#' \code{"gp"} and \code{sample_gp_hyperparameters} is \code{TRUE}.
#' @param auto_tune Logical; if \code{TRUE} (the default), every proposal
#' window above is adapted during the first \code{n_burn} iterations via
#' Robbins-Monro diminishing adaptation, instead of staying fixed at the
#' value passed in.
#' @param n_burn Number of iterations treated as burn-in for proposal-window
#' adaptation; ignored if \code{auto_tune} is \code{FALSE}. Defaults to half
#' of \code{n_iter} (resolved by \code{\link{runBatchMix}}, not here, since
#' \code{n_iter} is not known at \code{control} construction time).
#' @return An object of class \code{batchmix_control}: a named list holding
#' the arguments above, validated for type/length only (the interaction
#' between a window and the fitted model - e.g. whether it is suspiciously
#' large - is checked later, once \code{type} is known).
#' @export
#' @examples
#' # Defaults, equivalent to not passing `control` at all
#' batchmixControl()
#'
#' # Tighter cluster-mean proposal, auto-tuning switched off
#' batchmixControl(mu_proposal_window = 0.1, auto_tune = FALSE)
batchmixControl <- function(mu_proposal_window = 0.5**2,
                            cov_proposal_window = 0.002,
                            r_proposal_window = 0.1,
                            sigma_proposal_window = 0.01,
                            m_proposal_window = 0.3**2,
                            S_proposal_window = 0.01,
                            t_df_proposal_window = 0.015,
                            gamma_proposal_window = 0.1,
                            eta_proposal_window = 0.1,
                            gp_hyperparameter_proposal_window = 0.1,
                            auto_tune = TRUE,
                            n_burn = NULL) {
  control <- list(
    mu_proposal_window = mu_proposal_window,
    cov_proposal_window = cov_proposal_window,
    r_proposal_window = r_proposal_window,
    sigma_proposal_window = sigma_proposal_window,
    m_proposal_window = m_proposal_window,
    S_proposal_window = S_proposal_window,
    t_df_proposal_window = t_df_proposal_window,
    gamma_proposal_window = gamma_proposal_window,
    eta_proposal_window = eta_proposal_window,
    gp_hyperparameter_proposal_window = gp_hyperparameter_proposal_window,
    auto_tune = auto_tune,
    n_burn = n_burn
  )

  window_fields <- setdiff(names(control), c("auto_tune", "n_burn"))
  for (field in window_fields) {
    value <- control[[field]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
      stop(
        sprintf("batchmixControl(): `%s` must be a single numeric value.", field),
        call. = FALSE
      )
    }
  }

  if (!is.logical(auto_tune) || length(auto_tune) != 1L || is.na(auto_tune)) {
    stop("batchmixControl(): `auto_tune` must be a single logical value (TRUE/FALSE).", call. = FALSE)
  }

  if (!is.null(n_burn) && (!is.numeric(n_burn) || length(n_burn) != 1L || is.na(n_burn) || n_burn < 0)) {
    stop("batchmixControl(): `n_burn` must be NULL or a single non-negative number.", call. = FALSE)
  }

  structure(control, class = "batchmix_control")
}

#' @rdname batchmixControl
#' @param x An object of class \code{batchmix_control}.
#' @param ... Unused; present for S3 method consistency.
#' @export
print.batchmix_control <- function(x, ...) {
  cat("<batchmix control>\n")
  cat(sprintf(
    "  auto_tune = %s, n_burn = %s\n",
    x$auto_tune, if (is.null(x$n_burn)) "NULL (defaults to n_iter / 2)" else x$n_burn
  ))
  cat("  Proposal windows (starting values if auto_tune = TRUE, fixed otherwise):\n")
  window_fields <- setdiff(names(x), c("auto_tune", "n_burn"))
  for (field in window_fields) {
    cat(sprintf("    %-35s %s\n", field, format(x[[field]])))
  }
  invisible(x)
}

# Internal helper backing the individual-proposal-window/auto_tune/n_burn ->
# `control = batchmixControl()` deprecation shim, mirroring
# .resolveDeprecatedNIter()'s new-wins-over-old precedent: if `control` was
# explicitly supplied, it wins over any individually-supplied deprecated
# argument for the same field (with a warning naming the ignored ones); if
# `control` was left at its default, individually-supplied deprecated
# arguments are merged onto it instead (also with a warning).
# control_missing: missing(control) from the caller.
# control_value: the caller's `control` argument, already evaluated.
# deprecated: a named list of the caller's individually-supplied deprecated
# arguments (the caller builds this with `if (!missing(x)) list(x = x)`,
# `c()`-combined - see batchSemiSupervisedMixtureModel()/runBatchMix()/
# fitBatchMix() for the pattern).
# fn_name: the calling function's name, for the warning message.
.resolveControlArgs <- function(control_missing, control_value, deprecated, fn_name) {
  if (control_missing) {
    control <- batchmixControl()
  } else {
    if (!inherits(control_value, "batchmix_control")) {
      stop("`control` must be created with batchmixControl().", call. = FALSE)
    }
    control <- control_value
  }

  if (length(deprecated) > 0) {
    field_list <- paste(sprintf("`%s`", names(deprecated)), collapse = ", ")

    if (control_missing) {
      warning(
        sprintf(
          paste0(
            "Passing %s directly to %s() is deprecated; pass them inside ",
            "`control = batchmixControl(...)` instead. %s will continue to ",
            "work for this release but may be removed in a future one."
          ),
          field_list, fn_name, field_list
        ),
        call. = FALSE
      )
      control[names(deprecated)] <- deprecated
    } else {
      warning(
        sprintf(
          "Both `control` and %s were supplied to %s(); using `control` and ignoring %s.",
          field_list, fn_name, field_list
        ),
        call. = FALSE
      )
    }
  }

  control
}
