#!/usr/bin/Rscript
# Internal helper backing the `R` -> `n_iter` deprecation shim on the
# package's main entry points (runBatchMix, batchSemiSupervisedMixtureModel,
# continueChain, continueChains). `R` is deliberately captured via `...`
# rather than as a normal trailing `R = NULL` formal: a real formal shifts
# positional matching for every OTHER unnamed argument once `R` is supplied
# by name (e.g. `f(X, R = 1000, thin, batch_vec, type)` would then bind
# `thin`'s value to `n_iter`, not `thin`), which silently breaks exactly the
# legacy call style this shim exists to support. Passing `...` through
# leaves ordinary positional matching of the real formals untouched no
# matter how `R` is supplied.
# n_iter_missing: missing(n_iter) from the caller.
# n_iter_value: the caller's n_iter value, or NULL if missing.
# dots: the caller's list(...).
# fn_name: the calling function's name, for the "unused argument" error
# message (mimics base R's own wording).
.resolveDeprecatedNIter <- function(n_iter_missing, n_iter_value, dots, fn_name) {
  has_R <- "R" %in% names(dots)
  extra <- setdiff(names(dots), "R")
  if (length(extra) > 0) {
    stop(
      sprintf(
        "unused argument%s in %s(): %s",
        if (length(extra) > 1) "s" else "", fn_name,
        paste(sprintf("`%s`", extra), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!has_R) {
    if (n_iter_missing) {
      stop("argument \"n_iter\" is missing, with no default", call. = FALSE)
    }
    return(n_iter_value)
  }

  if (n_iter_missing) {
    warning(
      "The `R` argument is deprecated and has been renamed to `n_iter`. ",
      "`R` will continue to work for this release but may be removed in a ",
      "future one - please update calling code to use `n_iter` instead.",
      call. = FALSE
    )
    return(dots$R)
  }

  warning(
    "Both `n_iter` and the deprecated `R` argument were supplied; using ",
    "`n_iter` and ignoring `R`.",
    call. = FALSE
  )
  n_iter_value
}
