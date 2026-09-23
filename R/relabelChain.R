#!/usr/bin/Rscript
#' @title Relabel an MCMC chain to correct for label switching
#' @description Mixture models are only identified up to a permutation of
#' the component labels ("label switching", see Stephens (2000), "Dealing
#' with label switching in mixture models", JRSS-B 62(4), and Celeux, Hurn &
#' Robert (2000), JASA 95(451)): component 1 in one iteration may correspond
#' to component 2 in the next, with no bearing on convergence. Averaging (or
#' taking the median of) raw, cluster-indexed quantities across iterations -
#' \code{means}, \code{covariance}, \code{mean_sum}, \code{cov_comb},
#' \code{weights}, \code{t_df}, \code{gamma}, \code{w_batch}, \code{alloc} -
#' is only valid once every iteration's labels have been aligned to a common
#' reference. \code{\link{assessConvergence}}/\code{\link{rankNormalizedRhat}}
#' sidestep this by working on the permutation-invariant complete-data
#' log-likelihood/BIC trace instead (see their documentation) rather than
#' raw parameters, and \code{\link{minVI}}/\code{\link{VI.lb}} sidestep it by
#' working on the posterior similarity matrix; this function instead
#' actually resolves the permutation, for callers (chiefly
#' \code{\link{processMCMCChain}}) that need per-component parameter point
#' estimates, not just a partition estimate or a convergence diagnostic.
#'
#' Every iteration's hard allocation vector (\code{samples}) is matched to a
#' single reference labelling (the last retained iteration) via the
#' Equivalence Classes Representatives approach (Papastamoulis & Iliopoulos,
#' 2010, "An artificial allocations based solution to the label switching
#' problem in Bayesian analysis of mixtures of distributions", JCGS 19(2)):
#' for each iteration, the permutation of component labels that maximises
#' agreement with the reference is found by solving a linear sum assignment
#' problem (Hungarian algorithm, \code{\link[clue]{solve_LSAP}}) on the K x K
#' overlap/contingency table between that iteration's labels and the
#' reference's. Any component permanently anchored by at least one
#' semi-supervised (\code{fixed}) item is, by construction, never actually
#' moved by this (that item's label never changes sweep to sweep, so its
#' component trivially has near-total overlap with itself in the
#' reference) - only components with no fixed anchor (fully unsupervised
#' fits, or "extra" components beyond the labelled classes) can be, and are,
#' genuinely relabelled.
#' @param mcmc_output Output from \code{\link{batchSemiSupervisedMixtureModel}}
#' (or an entry of the list from \code{\link{runMCMCChains}}/
#' \code{\link{fitBatchMix}}), already burn-in trimmed if desired - every
#' array's third (iteration) dimension is relabelled in place, so a burn-in
#' should ordinarily be applied first (see \code{\link{processMCMCChain}}).
#' @param K_max The number of components; defaults to \code{mcmc_output$K_max}.
#' @return \code{mcmc_output} with \code{samples} and every present
#' cluster-indexed array (\code{means}, \code{covariance}, \code{mean_sum},
#' \code{cov_comb}, \code{weights}, \code{t_df}, \code{gamma}, \code{w_batch},
#' \code{alloc}) permuted, iteration by iteration, onto the common reference
#' labelling. Fields that are absent (e.g. \code{t_df} for a non-MVT fit) are
#' left untouched.
#' @export
relabelChain <- function(mcmc_output, K_max = mcmc_output[["K_max"]]) {
  # [[ ]], not $, throughout this function: mcmc_output can have dozens of
  # similarly-prefixed fields (e.g. `t_df_proposal_window`), and `$`'s
  # partial matching silently returns one of those instead of NULL when
  # the short name has no exact field - confirmed in practice: for a
  # non-MVT fit (no real `t_df` array), `mcmc_output$t_df` silently
  # returned `t_df_proposal_window`'s value instead of NULL, so
  # `has_t_df` below was wrongly TRUE and every downstream `t_df` access
  # operated on a scalar, not an (n_saved x K) matrix. `[[ ]]` only ever
  # matches the exact name, returning NULL when a field is genuinely
  # absent, exactly as intended.
  samples <- mcmc_output[["samples"]] # n_saved x N hard-label draws
  n_saved <- if (is.matrix(samples)) nrow(samples) else 0L

  # Nothing to relabel with fewer than 2 components, or fewer than 2 saved
  # iterations to align to a reference.
  if (n_saved < 2 || is.null(K_max) || K_max < 2) {
    return(mcmc_output)
  }

  # `samples` (and hence `alloc`'s allocation) holds label VALUES that are
  # 0-indexed (0..K_max-1) - the raw C++ sampler's own convention (labels
  # are an arma::uword there; see e.g. class_record.row(save_int) =
  # my_sampler.labels.t() in sampleSemisupervisedMVN.cpp et al.), unlike
  # every other cluster-indexed array here (means/covariance/.../alloc's
  # own K-axis), which are addressed by ordinary 1-indexed R array
  # POSITIONS, not stored values. Shift to 1-indexed only for this
  # function's own label-matching/permutation-application arithmetic
  # (which needs valid 1-indexed R vector indices throughout), and shift
  # back before returning - getting this wrong silently corrupts
  # `perm[samples[t, ]]` whenever a stored label is 0, since R indexing
  # with 0 drops that element rather than erroring loudly at the source
  # (it surfaces downstream instead, as a replacement-length mismatch).
  samples <- samples + 1L

  P <- mcmc_output[["P"]]
  B <- mcmc_output[["B"]]

  has_means <- !is.null(mcmc_output[["means"]])
  has_cov <- !is.null(mcmc_output[["covariance"]])
  has_weights <- !is.null(mcmc_output[["weights"]])
  has_t_df <- !is.null(mcmc_output[["t_df"]])
  has_mean_sum <- !is.null(mcmc_output[["mean_sum"]])
  has_cov_comb <- !is.null(mcmc_output[["cov_comb"]])
  has_gamma <- !is.null(mcmc_output[["gamma"]]) && isTRUE(mcmc_output[["include_interaction"]])
  has_w_batch <- !is.null(mcmc_output[["w_batch"]]) &&
    !is.null(mcmc_output[["weight_prior_type"]]) && mcmc_output[["weight_prior_type"]] > 0
  has_alloc <- !is.null(mcmc_output[["alloc"]])

  z_ref <- samples[n_saved, ]
  identity_perm <- seq_len(K_max)

  for (t in seq_len(n_saved)) {
    perm <- .hungarianComponentPermutation(samples[t, ], z_ref, K_max)

    # Already aligned (e.g. the reference iteration itself, or a component
    # arrangement that already agrees with it) - nothing to permute.
    if (all(perm == identity_perm)) {
      next
    }

    # inv_perm[j] is the ORIGINAL (this iteration's) component index that
    # should be written into aligned column/slice j, i.e. the inverse of
    # `perm` (perm[k] = reference component matched to this iteration's k).
    inv_perm <- integer(K_max)
    inv_perm[perm] <- identity_perm

    samples[t, ] <- perm[samples[t, ]]

    if (has_means) {
      mcmc_output[["means"]][, , t] <- mcmc_output[["means"]][, inv_perm, t]
    }

    if (has_cov) {
      # Contiguous P-wide blocks, one per cluster (block k = columns
      # ((k-1)*P+1):(k*P) - see cov_saved in sampleSemisupervisedMVN.cpp
      # et al.: a K-cube of P x P slices reshaped cluster-major).
      block_cols <- as.vector(vapply(
        inv_perm, function(k) ((k - 1) * P + 1):(k * P),
        integer(P)
      ))
      mcmc_output[["covariance"]][, , t] <- mcmc_output[["covariance"]][, block_cols, t]
    }

    if (has_weights) {
      mcmc_output[["weights"]][t, ] <- mcmc_output[["weights"]][t, inv_perm]
    }

    if (has_t_df) {
      mcmc_output[["t_df"]][t, ] <- mcmc_output[["t_df"]][t, inv_perm]
    }

    if (has_mean_sum) {
      # mean_sum's column for (cluster k, batch b) is (k-1)*B + b (the C++
      # sampler indexes both by kb = k * B + b) - contiguous B-wide blocks,
      # one per cluster.
      block_cols <- as.vector(vapply(
        inv_perm, function(k) ((k - 1) * B + 1):(k * B),
        integer(B)
      ))
      mcmc_output[["mean_sum"]][, , t] <- mcmc_output[["mean_sum"]][, block_cols, t]
    }

    if (has_cov_comb) {
      # cov_comb's block for (cluster k, batch b) is P columns wide at the
      # same kb = (k-1)*B + b block position as mean_sum, i.e. contiguous
      # (B*P)-wide blocks, one per cluster.
      block_cols <- as.vector(vapply(
        inv_perm, function(k) ((k - 1) * B * P + 1):(k * B * P),
        integer(B * P)
      ))
      mcmc_output[["cov_comb"]][, , t] <- mcmc_output[["cov_comb"]][, block_cols, t]
    }

    if (has_gamma) {
      # gamma is a genuine (P, K, B) cube flattened by Armadillo's own
      # column-major layout, giving column (batch b, cluster k) = (b-1)*K + k
      # - i.e. K-wide contiguous blocks, one per BATCH (not per cluster), so
      # each of the B blocks needs its K columns permuted independently.
      g <- mcmc_output[["gamma"]][, , t]
      for (b in seq_len(B)) {
        cols <- ((b - 1) * K_max + 1):(b * K_max)
        g[, cols] <- g[, cols, drop = FALSE][, inv_perm, drop = FALSE]
      }
      mcmc_output[["gamma"]][, , t] <- g
    }

    if (has_w_batch) {
      # w_batch is (B, K, n_saved) - cluster is the 2nd (column) axis.
      mcmc_output[["w_batch"]][, , t] <- mcmc_output[["w_batch"]][, inv_perm, t]
    }

    if (has_alloc) {
      # alloc is (N, K, n_saved) - cluster is the 2nd (column) axis.
      mcmc_output[["alloc"]][, , t] <- mcmc_output[["alloc"]][, inv_perm, t]
    }
  }

  mcmc_output[["samples"]] <- samples - 1L
  mcmc_output
}

# Equivalence Classes Representatives (Papastamoulis & Iliopoulos, 2010):
# find the permutation of z's component labels that maximises agreement
# with z_ref, via the Hungarian algorithm on the K x K overlap table.
# Returns perm such that perm[k] is the reference component matched to z's
# component k.
.hungarianComponentPermutation <- function(z, z_ref, K_max) {
  overlap <- matrix(0L, K_max, K_max)
  for (k in seq_len(K_max)) {
    idx_k <- which(z == k)
    if (length(idx_k) == 0) {
      next
    }
    overlap[k, ] <- tabulate(z_ref[idx_k], nbins = K_max)
  }
  as.integer(clue::solve_LSAP(overlap, maximum = TRUE))
}
