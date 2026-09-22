#!/usr/bin/Rscript
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
