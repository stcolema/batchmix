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
    "batchmix ", as.character(utils::packageVersion(pkgname)), ": this release ",
    "changes defaults and results, so existing analyses may differ.\n",
    "\n",
    "  Defaults (workflow-breaking):\n",
    "  * Partial pooling is now ON by default. Batches get their own mixture ",
    "weights (batch_weight_prior = \"partial_pooling\"; \"global\" for one batch) ",
    "and the concentration of the batch-scale prior is estimated ",
    "(sample_s_scale = TRUE). `alpha`/`concentration` are unused under partial ",
    "pooling, so an over-large K_max is no longer pushed towards empty ",
    "components by a sparse Dirichlet prior. For the old behaviour use ",
    "batch_weight_prior = \"global\", sample_s_scale = FALSE.\n",
    "\n",
    "  Results that changed:\n",
    "  * batch_corrected_data / inferred_dataset are now the posterior mean of ",
    "each item's batch-free signal (shrunk towards its cluster mean), the exact ",
    "inverse of the fitted batch model, not a rescaled residual.\n",
    "  * Proposal-window auto-tuning moved cov/S/sigma/t_df windows the wrong way ",
    "before; `mu_proposal_window` is now in units of each cluster's spread. ",
    "MVN_LKJ/MVN_MIXED previously targeted the wrong prior on the marginal ",
    "SDs. Fixed-seed output differs from earlier versions.\n",
    "  * Compare models with calcBICM() (best-draw BIC), not the per-iteration ",
    "BIC trace. Unsupervised point estimates now use minVI() on the posterior ",
    "similarity matrix (the `salso` dependency is gone; `nCores` is ignored).\n",
    "\n",
    "  Earlier changes: fitBatchMix() is the main entry point (runMCMCChains() is ",
    "a deprecated alias) and the iteration count `R` is now `n_iter`.\n",
    "\n",
    "  Please read the vignettes before re-running old analyses: ",
    "browseVignettes(\"batchmix\") - start with \"batchmix_workflow\", then ",
    "\"batch_weight_priors\" (what partial pooling changes) and ",
    "\"covariance_models\". See also news(package = \"batchmix\"). Silence this ",
    "message with suppressPackageStartupMessages(library(batchmix))."
  )
}
