#!/usr/bin/Rscript

# Column names referenced inside ggplot2::aes()/dplyr-style pipelines
# (tidyr::pivot_longer() output columns, or the .data pronoun itself) -
# these are resolved via non-standard evaluation at plot time, not as
# ordinary R variables, so R CMD check's static analysis cannot see where
# they come from and flags them as "no visible binding for global
# variable" without this declaration. This is the standard fix for
# NSE-heavy plotting code (see "Writing R Extensions" S1.6.4); it does not
# affect behaviour, only silences a check NOTE about code that is correct.
utils::globalVariables(c(
  ".data", "Acceptance_rate", "Chain", "Iteration", "Parameter", "iteration", "value"
))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "batchmix: recent changes to be aware of if you are updating from an ",
    "earlier version -\n",
    "  * The main entry point is now fitBatchMix() (fits multiple chains and ",
    "reports Rhat/ESS automatically); runMCMCChains() is now a deprecated ",
    "alias for it.\n",
    "  * The iteration-count argument `R` has been renamed to `n_iter` on ",
    "every entry point (runBatchMix(), batchSemiSupervisedMixtureModel(), ",
    "fitBatchMix(), continueChain(), continueChains()).\n",
    "  Both `runMCMCChains()` and the `R` argument still work for this ",
    "release (with a warning) but may be removed in a future one - please ",
    "update calling code. Use suppressPackageStartupMessages(library(batchmix)) ",
    "to silence this notice."
  )
}
