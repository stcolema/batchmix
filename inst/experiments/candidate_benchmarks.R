# Basic batchmix feasibility check on candidate benchmark datasets.
#
#   Rscript candidate_benchmarks.R <dataset> [n_iter] [out_dir]
#   dataset: gas | ifcb | tcga   (see the loaders below; each needs its data downloaded first)
#
# For every dataset the same report is produced: size (N, batches B, classes K, features P),
# class overlap (within-batch and pooled linear discriminant accuracy on the chosen features),
# whether the model fits (Rhat of the complete likelihood, acceptance rates, run time),
# classification of unlabelled items, and recovery of per-batch class proportions from the
# allocations. Models: global weights, partial pooling, partial pooling + interaction, and
# (when the dataset has a batch coordinate) a GP random walk on the weights.
suppressMessages({library(batchmix); library(MASS)})
args <- commandArgs(TRUE)
dataset <- if (length(args) >= 1) args[1] else "gas"
n_iter <- if (length(args) >= 2) as.integer(args[2]) else 2000L
out_dir <- if (length(args) >= 3) args[3] else "candidate_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
data_dir <- Sys.getenv("CANDIDATE_DATA_DIR", unset = "/tmp/claude-0/cand")
set.seed(42)

# ---- generic assessment -------------------------------------------------------------------------
lda_cv <- function(X, g, folds = 5) {
  g <- factor(g); ok <- table(g) >= folds; if (sum(ok) < 2) return(NA)
  keep <- g %in% names(ok)[ok]; X <- X[keep, , drop = FALSE]; g <- droplevels(g[keep])
  f <- sample(rep(seq_len(folds), length.out = nrow(X)))
  mean(sapply(seq_len(folds), function(k) {
    m <- tryCatch(lda(X[f != k, , drop = FALSE], g[f != k]), error = function(e) NULL)
    if (is.null(m)) return(NA)
    mean(predict(m, X[f == k, , drop = FALSE])$class == g[f == k])
  }), na.rm = TRUE)
}

assess <- function(name, X, batch, truth, coords = NULL, K, fixed_frac = 0.1, specs = NULL) {
  B <- length(unique(batch)); N <- nrow(X); P <- ncol(X)
  cat("\n", strrep("=", 78), "\n", name, ": N =", N, "| batches B =", B, "| classes K =", K, "| features P =", P, "\n", strrep("=", 78), "\n", sep = " ")
  cat("batch sizes (min/median/max):", range(table(batch))[1], median(table(batch)), range(table(batch))[2], "\n")
  cat("class counts:", table(truth), "\n")
  # class overlap
  within <- sapply(sort(unique(batch)), function(b) lda_cv(X[batch == b, , drop = FALSE], truth[batch == b]))
  pooled <- lda_cv(X, truth)
  cat(sprintf("overlap: within-batch LDA CV accuracy mean %.3f (range %.3f-%.3f); pooled across batches, no correction %.3f\n",
              mean(within, na.rm = TRUE), min(within, na.rm = TRUE), max(within, na.rm = TRUE), pooled))
  # between-batch drift in units of the within-class sd, on the first feature
  cm <- tapply(X[, 1], list(batch, truth), median); csd <- tapply(X[, 1], truth, sd)
  cat(sprintf("drift (feature 1): median over classes of the range of per-batch class medians / pooled class sd = %.2f\n",
              median(apply(cm, 2, function(z) diff(range(z, na.rm = TRUE))) / csd, na.rm = TRUE)))
  # supervision: fixed_frac of each class within each batch keeps its label
  fixed <- integer(N)
  for (b in unique(batch)) for (k in unique(truth)) {
    i <- which(batch == b & truth == k)
    if (length(i) > 0) fixed[sample(i, max(1, round(fixed_frac * length(i))))] <- 1L
  }
  cls <- sort(unique(truth)); lab_idx <- match(truth, cls) - 1L
  mu <- sapply(seq_along(cls), function(k) colMeans(X[fixed == 1 & lab_idx == k - 1L, , drop = FALSE]))
  S <- cov(X[fixed == 1, , drop = FALSE]) + diag(1e-6, P)
  init <- apply(X, 1, function(x) which.min(colSums(solve(S, x - mu) * (x - mu))) - 1L)
  init[fixed == 1] <- lab_idx[fixed == 1]
  unl <- fixed == 0
  true_prop <- table(factor(batch, levels = sort(unique(batch))), factor(lab_idx, levels = seq_len(K) - 1L)) /
    as.vector(table(batch))
  cat("fixed (labelled) items:", sum(fixed), "| unlabelled:", sum(unl), "| initial nearest-mean accuracy on unlabelled:",
      round(mean(init[unl] == lab_idx[unl]), 3), "\n")
  batch1 <- match(batch, sort(unique(batch)))
  if (is.null(specs)) specs <- list(global = list(batch_weight_prior = "global"),
                                    pp = list(batch_weight_prior = "partial_pooling"),
                                    pp_int = list(batch_weight_prior = "partial_pooling", include_interaction = TRUE))
  if (!is.null(coords)) specs$rw1 <- list(batch_weight_prior = "gp", gp_kernel = "rw1")
  rows <- list(); fits <- list()
  for (nm in names(specs)) {
    a <- c(list(X, 4, n_iter, max(1, n_iter %/% 200), batch1, "MVN", initial_labels = init, fixed = fixed,
                control = batchmixControl(n_burn = n_iter %/% 3), K_max = K), specs[[nm]])
    if (identical(specs[[nm]]$batch_weight_prior, "gp")) { a$batch_coordinates <- coords; a$sample_gp_hyperparameters <- TRUE }
    t0 <- Sys.time()
    o <- tryCatch(suppressMessages(do.call(fitBatchMix, a)), error = function(e) e)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(o, "error")) { rows[[nm]] <- data.frame(model = nm, error = conditionMessage(o)); next }
    best <- getBestChain(o); pc <- processMCMCChain(best, n_iter %/% 2)
    pred <- pc$pred - 1L
    est_prop <- table(factor(batch1, levels = seq_len(B)), factor(pred, levels = seq_len(K) - 1L)) / as.vector(table(batch1))
    # per-class recall on unlabelled items
    rec <- sapply(seq_len(K) - 1L, function(k) mean(pred[unl & lab_idx == k] == k))
    rows[[nm]] <- data.frame(model = nm, rhat = attr(o, "convergence")$rhat, secs = round(secs),
                             acc_unlabelled = mean(pred[unl] == lab_idx[unl]), min_class_recall = min(rec, na.rm = TRUE),
                             prop_MAE = mean(abs(est_prop - true_prop)),
                             acc_mu = mean(best$mu_acceptance_rate), acc_m = mean(best$m_acceptance_rate))
    fits[[nm]] <- list(pc = pc, rec = rec, est_prop = est_prop)
    cat(sprintf("  %-7s done in %.0fs\n", nm, secs))
  }
  res <- do.call(rbind, lapply(rows, function(r) r))
  cat("\nfit summary (best chain by BICM; Rhat > 1.01 means provisional):\n"); print(res, digits = 3, row.names = FALSE)
  if (length(fits)) { cat("\nper-class recall on unlabelled items, first model with a fit:\n"); print(round(fits[[1]]$rec, 3)) }
  list(res = res, within = within, pooled = pooled, true_prop = true_prop, fits = lapply(fits, function(f) f[c("rec", "est_prop")]))
}

# ---- loaders ----------------------------------------------------------------------------------------
load_gas <- function(max_per_cell = 200, conc_range = NULL) {
  # UCI gas sensor array drift (Vergara et al. 2012): 16 sensors x 8 features, 6 gases, 10 batches
  d <- readRDS(file.path(data_dir, "gas.rds"))
  sl <- function(x) sign(x) * log1p(abs(x))               # signed log: 0.2% of the first feature is <= 0
  if (!is.null(conc_range)) d <- d[d$conc > conc_range[1] & d$conc <= conc_range[2], ]   # one concentration band
  X <- sapply(1:4, function(s) sl(d[[paste0("X", (s - 1) * 8 + 1)]]))   # steady-state response, sensors 1-4
  colnames(X) <- paste0("sensor", 1:4)
  idx <- unlist(lapply(split(seq_len(nrow(d)), list(d$batch, d$gas)), function(i) if (length(i) > max_per_cell) sample(i, max_per_cell) else i))
  idx <- sort(idx)
  # approximate month of each batch (UCI description; recalled, not checked against the files)
  months <- c(1.5, 6.5, 12, 14.5, 16, 18.5, 21, 22.5, 27, 36)
  list(X = X[idx, ], batch = d$batch[idx], truth = d$gas[idx], coords = months, K = 6)
}

load_tcga <- function(k = 2, min_plate = 15) {
  # TCGA-BRCA: RSEM mRNA (cBioPortal) + IHC/FISH receptor status; plate from the GDC aliquot barcode
  suppressMessages(library(jsonlite))
  td <- file.path(data_dir, "tcga")
  e <- fromJSON(file.path(td, "expr.json")); sym <- c(`2099` = "ESR1", `5241` = "PGR", `2064` = "ERBB2", `2886` = "GRB7",
    `2625` = "GATA3", `3169` = "FOXA1", `4288` = "MKI67", `3852` = "KRT5", `3872` = "KRT17")
  e$gene <- sym[as.character(e$entrezGeneId)]
  W <- reshape(e[, c("sampleId", "gene", "value")], idvar = "sampleId", timevar = "gene", direction = "wide")
  names(W) <- sub("value\\.", "", names(W))
  cl <- fromJSON(file.path(td, "clin_pat.json"))
  cw <- reshape(cl[cl$clinicalAttributeId %in% c("ER_STATUS_BY_IHC", "PR_STATUS_BY_IHC", "IHC_HER2", "HER2_FISH_STATUS"),
                   c("patientId", "clinicalAttributeId", "value")], idvar = "patientId", timevar = "clinicalAttributeId", direction = "wide")
  names(cw) <- sub("value\\.", "", names(cw))
  al <- read.csv(file.path(td, "aliquots.csv")); al <- al[al$sample_type == "Primary Tumor", ]
  al$sampleId <- sub("[A-Z]$", "", al$sample); al$plate <- sapply(strsplit(al$aliquot, "-"), `[`, 6)
  al <- al[!duplicated(al$sampleId), c("sampleId", "plate", "patient")]
  d <- merge(merge(W, al, by = "sampleId"), cw, by.x = "patient", by.y = "patientId")
  # HER2: positive if IHC or FISH is positive; negative if at least one is negative and none positive
  # (equivocal / indeterminate / blank results are treated as absent)
  hp <- (d$IHC_HER2 %in% "Positive") | (d$HER2_FISH_STATUS %in% "Positive")
  hn <- (d$IHC_HER2 %in% "Negative") | (d$HER2_FISH_STATUS %in% "Negative")
  her2 <- ifelse(hp, "pos", ifelse(hn, "neg", NA))
  er <- d$ER_STATUS_BY_IHC; pr <- d$PR_STATUS_BY_IHC
  if (k == 2) {
    cls <- ifelse(er == "Positive", 1L, ifelse(er == "Negative", 0L, NA))
    feats <- c("ESR1", "PGR", "GATA3", "FOXA1", "ERBB2", "GRB7")
  } else {
    cls <- ifelse(her2 == "pos", 2L, ifelse(her2 == "neg" & (er == "Positive" | pr == "Positive"), 1L,
           ifelse(her2 == "neg" & er == "Negative" & pr == "Negative", 0L, NA)))   # 0 TNBC, 1 HR+/HER2-, 2 HER2+
    feats <- c("ESR1", "PGR", "ERBB2", "GRB7", "KRT5", "KRT17")
  }
  X <- log2(as.matrix(d[, feats]) + 1)
  keep <- !is.na(cls) & complete.cases(X)
  d <- d[keep, ]; X <- X[keep, ]; cls <- cls[keep]
  big <- names(which(table(d$plate) >= min_plate)); keep <- d$plate %in% big
  d <- d[keep, ]; X <- X[keep, ]; cls <- cls[keep]
  plates <- sort(unique(d$plate))
  # plate codes are issued sequentially, so their rank is used as a time-like coordinate (assumption, not checked)
  list(X = X, batch = match(d$plate, plates), truth = cls, coords = seq_along(plates), K = k)
}

load_ifcb <- function(K = 5, npc = 5) {
  # WHOI/MVCO IFCB training samples (2006-2008, one instrument): 512 deep features per image, manual class label.
  files <- list.files(file.path(data_dir, "ifcb", "train_sub"), full.names = TRUE)
  L <- lapply(files, function(f) { z <- read.csv(gzfile(f)); z$sample <- sub("\\.csv\\.gz$", "", basename(f)); z })
  d <- do.call(rbind, L)
  top <- as.integer(names(sort(table(d$class), decreasing = TRUE))[seq_len(K)])
  d <- d[d$class %in% top, ]
  F <- as.matrix(d[, as.character(0:511)])
  pc <- prcomp(F, center = TRUE, scale. = FALSE); X <- pc$x[, seq_len(npc)]
  yd <- as.integer(sub(".*_(\\d{4})_(\\d{3})_.*", "\\2", d$sample)); yr <- as.integer(sub(".*_(\\d{4})_.*", "\\1", d$sample))
  mo <- (yr - 2006) * 12 + ceiling(yd / 30.5)                      # month index since January 2006
  list(X = X, batch = match(mo, sort(unique(mo))), truth = match(d$class, top), coords = sort(unique(mo)), K = K,
       explained = summary(pc)$importance[3, npc])
}

fixed_frac <- as.numeric(Sys.getenv("FIXED_FRAC", unset = "0.1"))
if (dataset == "gas_band") {
  g <- load_gas(conc_range = c(25, 100))
  out <- assess("UCI gas drift, concentration band (25, 100] ppmv, sensors 1-4", g$X, g$batch, g$truth, g$coords, g$K, fixed_frac = fixed_frac)
}
if (dataset == "gas") {
  g <- load_gas()
  out <- assess("UCI gas sensor array drift (sensors 1-4, steady-state response)", g$X, g$batch, g$truth, g$coords, g$K, fixed_frac = fixed_frac)
}
if (dataset %in% c("tcga2", "tcga3")) {
  g <- load_tcga(if (dataset == "tcga2") 2 else 3)
  out <- assess(paste("TCGA-BRCA mRNA,", if (dataset == "tcga2") "ER status (K=2)" else "TNBC / HR+ HER2- / HER2+ (K=3)"),
                g$X, g$batch, g$truth, g$coords, g$K)
}
if (dataset == "ifcb") {
  g <- load_ifcb()
  cat("PCA: variance explained by", ncol(g$X), "components:", round(g$explained, 3), "\n")
  out <- assess("IFCB plankton (top 5 taxa, 5 PCs of the 512 CNN features, batch = month)", g$X, g$batch, g$truth, g$coords, g$K)
}
saveRDS(out, file.path(out_dir, paste0(dataset, "_results.rds")))
