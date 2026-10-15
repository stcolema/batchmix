# Parallel execution of chains.
#
# Chains are independent, so they can run on separate cores. The package never
# starts workers unless asked (`n_cores` or `options(batchmix.cores)`; the default
# is one core), and, as `data.table` does, never uses more than two cores while
# a package is being checked (`_R_CHECK_LIMIT_CORES_`, which `R CMD check
# --as-cran` sets), as the CRAN policy requires.

# Internal: the number of cores to use for n_chains chains
.batchmixResolveCores <- function(n_cores, n_chains) {
  if (is.null(n_cores)) {
    n_cores <- getOption("batchmix.cores", 1L)
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1 || is.na(n_cores) || n_cores < 1) {
    stop("`n_cores` must be a single positive integer.", call. = FALSE)
  }
  n_cores <- as.integer(floor(n_cores))
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))) {
    n_cores <- min(n_cores, 2L)
  }
  available <- parallel::detectCores(logical = TRUE)
  if (!is.na(available)) {
    n_cores <- min(n_cores, as.integer(available))
  }
  max(1L, min(n_cores, as.integer(n_chains)))
}

# Internal: one L'Ecuyer-CMRG stream per chain, so that chains are reproducible
# (and independent) whatever the number of cores. The seed is drawn from the
# user's generator, so set.seed() before the call fixes the result, and the
# user's generator kind and state are restored on exit.
.batchmixStreamSeeds <- function(n) {
  seed <- sample.int(.Machine$integer.max, 1L)
  envir <- globalenv()
  old_seed <- get(".Random.seed", envir = envir)
  old_kind <- RNGkind()
  on.exit({
    suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
    assign(".Random.seed", old_seed, envir = envir)
  }, add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  seeds <- vector("list", n)
  seeds[[1]] <- get(".Random.seed", envir = envir)
  for (i in seq_len(n)[-1]) {
    seeds[[i]] <- parallel::nextRNGStream(seeds[[i - 1]])
  }
  seeds
}

# Internal: run_chain(ii) for each chain ii, as a list. On one core this is
# lapply() and draws from the user's generator exactly as before. On several, each
# chain gets its own stream (see .batchmixStreamSeeds()), so the result follows
# set.seed() and does not depend on the number of cores, nor on whether workers are
# forked (Unix-alikes) or sockets (Windows); `options(batchmix.parallel_type)` can
# force either. It differs from the serial result for the same seed.
.batchmixRunChains <- function(n_chains, n_cores, run_chain) {
  if (n_cores <= 1) {
    return(lapply(seq_len(n_chains), run_chain))
  }

  seeds <- .batchmixStreamSeeds(n_chains)
  worker <- function(ii) {
    assign(".Random.seed", seeds[[ii]], envir = globalenv())
    run_chain(ii)
  }

  type <- getOption("batchmix.parallel_type", if (.Platform$OS.type == "windows") "PSOCK" else "FORK")
  if (identical(type, "FORK")) {
    # A failed chain is reported below; mclapply() also warns about it
    out <- withCallingHandlers(
      parallel::mclapply(
        seq_len(n_chains), worker,
        mc.cores = n_cores, mc.preschedule = FALSE, mc.set.seed = FALSE
      ),
      warning = function(w) {
        if (grepl("encountered errors in user code|did not deliver|resulted in an error", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  } else {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Workers must find this package wherever the session does
    parallel::clusterCall(cl, function(lib_paths) .libPaths(lib_paths), .libPaths())
    out <- parallel::parLapply(cl, seq_len(n_chains), worker)
  }

  failed <- which(vapply(out, function(x) inherits(x, "try-error") || is.null(x), logical(1)))
  if (length(failed) > 0) {
    reasons <- vapply(out[failed], function(x) {
      if (inherits(x, "try-error")) trimws(as.character(x)) else "no result returned"
    }, character(1))
    stop(
      "Chain", if (length(failed) > 1) "s " else " ", paste(failed, collapse = ", "),
      " failed: ", reasons[1],
      call. = FALSE
    )
  }
  out
}
