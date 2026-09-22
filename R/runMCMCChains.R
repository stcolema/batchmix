#!/usr/bin/Rscript
#' @title Run MCMC Chains (deprecated)
#' @description Deprecated alias for \code{\link{fitBatchMix}} - kept for
#' backwards compatibility only, and will be removed in a future release.
#' Every argument is forwarded unchanged (including the deprecated \code{R}
#' argument, which \code{fitBatchMix()} itself still handles); see
#' \code{?fitBatchMix} for the full, current documentation.
#' @param ... Forwarded to \code{\link{fitBatchMix}}.
#' @returns See \code{\link{fitBatchMix}}.
#' @export
runMCMCChains <- function(...) {
  warning(
    "runMCMCChains() is deprecated and has been renamed to fitBatchMix(); ",
    "please update calling code - runMCMCChains() will continue to work for ",
    "this release but may be removed in a future one.",
    call. = FALSE
  )
  fitBatchMix(...)
}
