# Bayesian workflow for the Stockholm SARS-CoV-2 serology data
# (De Castro Dopico et al., J Intern Med 2021; chr1swallace/seroprevalence-paper),
# with a batch x cluster interaction term and checking by posterior predictive
# checks (PPCs), held-out controls and the one repeatedly measured individual
# ("Patient 4").
#
#   Rscript stockholm_workflow.R /path/to/seroprevalence-paper [n_iter] [out_dir]
#
# Steps (after Gelman et al., 2020, "Bayesian Workflow", arXiv:2011.01808):
#   0  data and design
#   1  prior predictive check
#   2  fit a small set of models (computation and convergence diagnostics)
#   3  posterior predictive checks (batch-level statistics, tails)
#   4  held-out controls (sensitivity, specificity, calibration)
#   5  Patient 4: out-of-sample check of the plate-to-plate technical variation
#   6  comparison and sensitivity to the weight prior
#   7  weekly prevalence from the allocations (benchmark: published estimates)
#
# Design choices (not given by the data):
#  * batch = week (21 batches). Weeks nest inside four assay runs, and plate
#    IDs are not in the file, so run and week effects are confounded.
#  * features: log10 of the raw OD (spike, RBD), not the plate-adjusted columns.
#  * controls have no week: historical negatives are assigned at random to the
#    weeks of the run they were measured in (run 201007 has no donor weeks and is
#    dropped); COVID-positives to any week. Half of the controls keep their label
#    (fixed), half are unlabelled and used for scoring. Controls enter each
#    batch's class counts, which biases the weight prior toward their mix, so the
#    prevalence reported is the fraction of DONOR samples allocated to the
#    positive class, not the batch weight.
#  * Patient 4 (24 replicate measurements of one person, 22 in the donor block and
#    2 in the pregnant block) is held out of every fit. The file has no plate IDs.
suppressMessages({library(batchmix); library(data.table)})
args <- commandArgs(TRUE)
path <- if (length(args) >= 1) args[1] else "seroprevalence-paper"
n_iter <- if (length(args) >= 2) as.integer(args[2]) else 3000L
out_dir <- if (length(args) >= 3) args[3] else "stockholm_workflow"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(2024)
thin <- max(1L, n_iter %/% 200L)
n_burn <- n_iter %/% 3L
pdf(file.path(out_dir, "workflow_figures.pdf"), width = 9, height = 6)
section <- function(s) cat("\n", strrep("=", 78), "\n", s, "\n", strrep("=", 78), "\n", sep = "")
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- 0. Data and design -----------------------------------------------------
section("0. Data and design")
e <- new.env(); load(file.path(path, "adjusted-data.RData"), e); m <- as.data.table(e$m)
est <- fread(file.path(path, "estimates.csv"))
p4 <- m[type == "Patient 4"]                       # before any de-duplication
m[, week := suppressWarnings(as.integer(ifelse(grepl("Wk", Sample.ID), sub(".*Wk([0-9]+).*", "\\1", Sample.ID), NA)))]
m <- m[type != "Patient 4"][!duplicated(paste(group, type, Sample.ID))]
run_of_week <- function(w) ifelse(w <= 21, "early", ifelse(w <= 25, "200702", ifelse(w <= 34, "200923", "201216")))
donors <- m[type == "Blood donors"]
weeks <- sort(unique(donors$week))
runs <- split(weeks, run_of_week(weeks))
neg <- m[type == "Historical controls"]
neg[, run := ifelse(group %in% c("Blood Donors", "Pregnant volunteers"), "early", group)]
neg <- neg[run %in% names(runs)]
neg[, week := vapply(run, function(r) sample(rep(runs[[r]], 2), 1), 1L)]
pos <- m[type == "COVID"]
pos[, week := sample(weeks, .N, replace = TRUE)]
pos[, truth := 1L]; neg[, truth := 0L]; donors[, truth := NA_integer_]
d <- rbindlist(list(donors, neg, pos), fill = TRUE)
d[, role := ifelse(is.na(truth), "donor", ifelse(runif(.N) < 0.5, "fixed", "holdout"))]
tr <- function(x) log10(pmax(x, 0.01))
X <- cbind(spike = tr(d$OD.spike), rbd = tr(d$OD.rbd))
X_p4 <- cbind(spike = tr(p4$OD.spike), rbd = tr(p4$OD.rbd))
batch <- match(d$week, weeks)
fixed <- as.integer(d$role == "fixed")
is_donor <- d$role == "donor"; is_hold <- d$role == "holdout"
init <- ifelse(fixed == 1, d$truth, as.integer(X[, 1] > log10(0.25)))
cat("N =", nrow(d), "| donors", sum(is_donor), "| fixed controls", sum(fixed), "| held-out controls", sum(is_hold),
    "| Patient 4 replicates (held out)", nrow(p4), "\n")

# ---- 1. Prior predictive check ----------------------------------------------
section("1. Prior predictive check (partial pooling, default priors)")
pri <- simulatePriorPredictive(X, batch - 1L, K = 2, type = "MVN", batch_weight_prior = "partial_pooling", n_datasets = 40)
obs_rng <- apply(X, 2, quantile, c(0.01, 0.99))
cat("observed 1% and 99% quantiles of log10 OD:\n"); print(round(obs_rng, 2))
sim_rng <- sapply(pri, function(s) as.vector(apply(s$X, 2, quantile, c(0.01, 0.99))))
cat("prior predictive: median [5%, 95%] across datasets of the same quantiles (rows: low/high x feature)\n")
print(round(t(apply(sim_rng, 1, quantile, c(0.5, 0.05, 0.95))), 2))
cat("share of prior datasets with a log10 OD beyond +-6 (a 10^6-fold OD):",
    mean(sapply(pri, function(s) any(abs(s$X) > 6))), "\n")
par(mfrow = c(1, 2))
for (j in 1:2) {
  plot(density(X[, j]), lwd = 2, main = paste("Prior predictive:", colnames(X)[j]), xlab = "log10 OD",
       xlim = range(c(X[, j], unlist(lapply(pri[1:20], function(s) quantile(s$X[, j], c(.01, .99)))))))
  for (s in pri[1:20]) lines(density(s$X[, j]), col = rgb(0, 0, 1, 0.25))
}

# ---- 2. Fit ---------------------------------------------------------------------
section("2. Models and computation")
specs <- list(
  pp         = list(batch_weight_prior = "partial_pooling", include_interaction = FALSE),
  pp_int     = list(batch_weight_prior = "partial_pooling", include_interaction = TRUE),
  rw1_int    = list(batch_weight_prior = "gp", include_interaction = TRUE, gp_kernel = "rw1"),
  global_int = list(batch_weight_prior = "global", include_interaction = TRUE),
  pp_int_sd25 = list(batch_weight_prior = "partial_pooling", include_interaction = TRUE, pp_mu_prior_sd = 2.5)
)
fits <- list()
for (nm in names(specs)) {
  s <- specs[[nm]]
  a <- c(list(X, 4, n_iter, thin, batch, "MVN", initial_labels = init, fixed = fixed,
              control = batchmixControl(n_burn = n_burn), K_max = 2), s)
  if (!is.null(s$batch_weight_prior) && s$batch_weight_prior == "gp") {
    a$batch_coordinates <- weeks; a$sample_gp_hyperparameters <- TRUE
  }
  t0 <- Sys.time()
  o <- suppressMessages(do.call(fitBatchMix, a))
  best <- getBestChain(o)
  pc <- processMCMCChain(best, n_iter %/% 2)
  conv <- attr(o, "convergence")
  fits[[nm]] <- list(o = o, best = best, pc = pc, rhat = conv$rhat, ess = conv$ess_bulk %||% NA,
                     secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  cat(sprintf("%-12s Rhat(complete likelihood) %.3f | secs %.0f | acceptance: mu %.2f cov %.2f m %.2f S %.2f%s\n",
              nm, conv$rhat, fits[[nm]]$secs, mean(best$mu_acceptance_rate), mean(best$cov_acceptance_rate),
              mean(best$m_acceptance_rate), mean(best$S_acceptance_rate),
              if (isTRUE(s$include_interaction)) sprintf(" gamma %.2f", mean(best$gamma_acceptance_rate)) else ""))
  if (isTRUE(s$include_interaction)) {
    cat("             interaction variance tau2_gamma per feature, posterior mean:",
        round(rowMeans(pc$tau2_interaction %||% matrix(NA, 2, 1)), 4), "\n")
  }
}
cat("\nRhat above 1.01 means these results are provisional; rerun with a larger n_iter.\n")
# The likelihood sees only mu_k + m_b, so a common shift of all batches can be traded against the class means; the
# prior on m (centred at 0) is the only thing that pins it. Check: Rhat of the positive-class mean as sampled and
# after adding back the mean shift (mu + mean_b m_b, which the likelihood does identify), and how much of the
# fitted shifts' size is a common offset rather than between-batch spread.
ridge <- do.call(rbind, lapply(names(fits), function(nm) {
  per <- lapply(fits[[nm]]$o, function(ch) {
    px <- processMCMCChain(ch, n_iter %/% 2); n <- dim(px$means)[3]
    pk <- apply(px$means[1, , , drop = FALSE], 3, which.max)
    mu <- sapply(seq_len(n), function(r) px$means[1, pk[r], r])
    cbind(raw = mu, centred = mu + colMeans(px$batch_shift[1, , ]),
          common = colMeans(px$batch_shift[1, , ]), between_sd = apply(px$batch_shift[1, , ], 2, sd))
  })
  g <- function(col) do.call(cbind, lapply(per, function(z) z[, col]))
  data.frame(model = nm, rhat_mu_pos_raw = rankNormalizedRhat(g("raw"))$rhat,
             rhat_mu_pos_plus_mean_shift = rankNormalizedRhat(g("centred"))$rhat,
             rms_common_shift = sqrt(mean(g("common")^2)), mean_between_batch_sd = mean(g("between_sd")))
}))
cat("\nIdentifiability of the common shift (spike feature):\n"); print(ridge, digits = 3)

# ---- 3. Posterior predictive checks ------------------------------------------
section("3. Posterior predictive checks (batch-level statistics)")
n_rep <- 40
thr <- quantile(X[d$role %in% c("fixed", "holdout") & d$truth == 0, 1], 0.99)   # 99th pct of control negatives, spike
stat_fun <- function(Xm, idx) {
  # per-week mean and sd of each feature, and per-week share above the negative-control threshold
  w <- d$week[idx]
  out <- list()
  for (j in 1:2) {
    out[[paste0("mean_", j)]] <- tapply(Xm[idx, j], w, mean)
    out[[paste0("sd_", j)]]   <- tapply(Xm[idx, j], w, sd)
  }
  out$above <- tapply(Xm[idx, 1] > thr, w, mean)
  out$q99_1 <- quantile(Xm[idx, 1], 0.99); out$q99_2 <- quantile(Xm[idx, 2], 0.99)
  out
}
obs_stat <- stat_fun(X, which(is_donor))
ppc_tab <- list()
for (nm in names(fits)) {
  reps <- simulatePosteriorPredictive(fits[[nm]]$best, batch - 1L, n_draws = n_rep, burn = n_iter %/% 2, seed = 1)
  rs <- lapply(reps, function(r) stat_fun(r$X, which(is_donor)))
  pv <- lapply(names(obs_stat), function(k) {
    o_k <- obs_stat[[k]]
    r_k <- sapply(rs, function(z) z[[k]])
    if (is.null(dim(r_k))) r_k <- matrix(r_k, nrow = 1)
    rowMeans(r_k >= as.numeric(o_k))        # P(T_rep >= T_obs) for each week (or the single statistic)
  })
  names(pv) <- names(obs_stat)
  ext <- sapply(pv, function(p) mean(p < 0.025 | p > 0.975))
  ppc_tab[[nm]] <- ext
  fits[[nm]]$pv <- pv
}
cat("Share of posterior predictive p-values outside [0.025, 0.975] (about 0.05 is expected if the model fits;\n",
    "p-values from one chain's draws, conditional on its sampled allocations):\n", sep = "")
print(round(do.call(rbind, lapply(ppc_tab, function(z) z[!names(z) %in% c("q99_1", "q99_2")])), 2))
cat("\nUpper-tail checks, P(T_rep >= T_obs) for the 99th percentile of donor log10 OD (extreme if < 0.025 or > 0.975):\n")
print(round(t(sapply(fits, function(f) c(spike = f$pv$q99_1, rbd = f$pv$q99_2))), 3))
par(mfrow = c(2, 2))
for (nm in c("pp_int", "rw1_int")) for (k in c("mean_1", "sd_1")) {
  plot(weeks, fits[[nm]]$pv[[k]], ylim = c(0, 1), pch = 19, xlab = "week", ylab = "P(T_rep >= T_obs)",
       main = paste(nm, k)); abline(h = c(0.025, 0.5, 0.975), lty = c(2, 1, 2))
}

# ---- 4. Held-out controls -------------------------------------------------------
section("4. Held-out controls: classification and calibration")
hold_tab <- do.call(rbind, lapply(names(fits), function(nm) {
  pc <- fits[[nm]]$pc
  pp <- pc$allocation_probability[, 2]
  # the positive class is the one with the higher spike mean in the posterior mean
  if (mean(pc$mean_est[1, 2]) < mean(pc$mean_est[1, 1])) pp <- pc$allocation_probability[, 1]
  y <- d$truth[is_hold]; p <- pmin(pmax(pp[is_hold], 1e-6), 1 - 1e-6)
  data.frame(model = nm, sens = mean(p[y == 1] > 0.5), spec = mean(p[y == 0] <= 0.5),
             brier = mean((p - y)^2), log_score = mean(y * log(p) + (1 - y) * log(1 - p)))
}))
print(hold_tab, digits = 3)

# ---- 5. Patient 4 -------------------------------------------------------------------
section("5. Patient 4: plate-to-plate technical variation, out of sample")
cat("Patient 4 is one person measured", nrow(p4), "times, so the replicate spread is technical (plate, well,\n",
    "run), not biological. The model's technical component for an individual in class k and batch b is\n",
    "m_b + e, e ~ N(0, diag((S_b - 1) * diag(Sigma_k))). Replicates of one person on a new plate each:\n",
    "simulate", nrow(p4), "new batches per posterior draw and compare the sample sd with the observed sd.\n", sep = " ")
obs_sd <- apply(X_p4, 2, sd)
cat("observed sd of log10 OD across the replicates (spike, rbd):", round(obs_sd, 3), "\n")
cat("for scale: sd of log10 OD among historical-negative controls:",
    round(apply(X[d$truth %in% 0, ], 2, sd), 3), "\n")
p4_tab <- list(); p4_draws <- list()
for (nm in names(fits)) {
  pc <- fits[[nm]]$pc
  P <- pc$P; n_saved <- dim(pc$means)[3]
  pos_k <- apply(pc$means[1, , , drop = FALSE], 3, which.max)
  dG <- sapply(seq_len(n_saved), function(r) diag(pc$covariance[, ((pos_k[r] - 1) * P + 1):(pos_k[r] * P), r, drop = TRUE]))
  Y <- Ysh <- Ye <- array(NA_real_, c(nrow(p4), P, n_saved))
  for (j in seq_len(nrow(p4))) {
    sd_j <- batchmix:::.predictiveShiftScaleDraws(pc, X)
    Ysh[j, , ] <- sd_j$shift
    Ye[j, , ] <- matrix(rnorm(P * n_saved), P) * sqrt(dG * (sd_j$scale - 1))
    Y[j, , ] <- Ysh[j, , ] + Ye[j, , ]
  }
  rep_sd <- apply(Y, c(2, 3), sd)                  # P x n_saved
  sh_sd <- apply(Ysh, c(2, 3), sd); e_sd <- apply(Ye, c(2, 3), sd)
  # Diagnostic: the same prediction with the shift spread taken from the fitted BETWEEN-batch sd of the shifts
  # (what a centred hierarchy would learn) instead of lambda_2, which is inflated by the common offset.
  between <- apply(pc$batch_shift, c(1, 3), sd)    # P x n_saved
  Yc <- array(NA_real_, c(nrow(p4), P, n_saved))
  for (j in seq_len(nrow(p4))) {
    sd_j <- batchmix:::.predictiveShiftScaleDraws(pc, X)
    Yc[j, , ] <- matrix(rnorm(P * n_saved), P) * between + matrix(rnorm(P * n_saved), P) * sqrt(dG * (sd_j$scale - 1))
  }
  rep_sd_c <- apply(Yc, c(2, 3), sd)
  p4_draws[[nm]] <- rep_sd
  p4_tab[[nm]] <- c(
    pred_sd_spike = median(rep_sd[1, ]), pred_sd_rbd = median(rep_sd[2, ]),
    p_spike = mean(rep_sd[1, ] >= obs_sd[1]), p_rbd = mean(rep_sd[2, ] >= obs_sd[2]),
    shift_part_spike = median(sh_sd[1, ]), noise_part_spike = median(e_sd[1, ]),
    centred_pred_sd_spike = median(rep_sd_c[1, ]), centred_pred_sd_rbd = median(rep_sd_c[2, ]),
    centred_p_spike = mean(rep_sd_c[1, ] >= obs_sd[1]), centred_p_rbd = mean(rep_sd_c[2, ] >= obs_sd[2]),
    fitted_week_sd_spike = sd(pc$shift_est[1, ]), fitted_week_sd_rbd = sd(pc$shift_est[2, ])
  )
}
cat("\npredicted sd across", nrow(p4), "plates (median over draws), P(sd_rep >= sd_obs), and the sd of the fitted week shifts:\n")
print(round(do.call(rbind, p4_tab), 3))
cat("'centred_*' replaces the shift spread by the fitted between-batch sd of the shifts (a diagnostic, not the model as fitted).\n")
cat("A p near 0 means the model predicts less plate-to-plate spread than Patient 4 shows (the model is\n",
    "too optimistic about batch noise); near 1, more than shown. The fitted week shifts average over the\n",
    "plates in a week, so they understate plate-level spread.\n", sep = "")
par(mfrow = c(1, 2))
for (j in 1:2) {
  hist(p4_draws[["pp_int"]][j, ], breaks = 40, col = "grey85", border = NA, main = paste("Patient 4 sd,", colnames(X)[j]),
       xlab = "sd of log10 OD across plates"); abline(v = obs_sd[j], lwd = 2, col = 2)
}

# ---- 6. Comparison and sensitivity ---------------------------------------------------
section("6. Comparison: best-draw BIC (BICM), PPC and held-out score")
cmp <- do.call(rbind, lapply(names(fits), function(nm) data.frame(
  model = nm, BICM = calcBICM(fits[[nm]]$best, n_iter %/% 2), rhat = fits[[nm]]$rhat,
  ppc_extreme = mean(unlist(ppc_tab[[nm]])), log_score_holdout = hold_tab$log_score[hold_tab$model == nm])))
print(cmp, digits = 4)
cat("BICM differences between weight priors are only meaningful relative to their parameter counts\n",
    "(partial pooling and the GP add structural parameters); the PPC and held-out columns carry no such penalty.\n", sep = "")

# ---- 7. Weekly prevalence ---------------------------------------------------------------
section("7. Weekly prevalence among donors (fraction allocated to the positive class)")
ref <- est[type == "Blood donors" & grepl("^Wk", Week)]; ref[, week := as.integer(sub("Wk", "", Week))]
ref <- ref[match(weeks, week)]
prev <- do.call(cbind, lapply(names(fits), function(nm) {
  pc <- fits[[nm]]$pc
  pos_k <- apply(pc$means[1, , , drop = FALSE], 3, which.max) - 1L     # 0-indexed positive class per draw
  smp <- pc$samples[, is_donor, drop = FALSE]
  frac <- sapply(seq_len(nrow(smp)), function(r) tapply(smp[r, ] == pos_k[r], d$week[is_donor], mean))
  setNames(data.frame(rowMeans(frac), apply(frac, 1, quantile, 0.025), apply(frac, 1, quantile, 0.975)),
           paste0(nm, c("", "_lo", "_hi")))
}))
tab7 <- cbind(week = weeks, published = ref$estimate, pub_lo = ref$lower.ci, pub_hi = ref$upper.ci, prev)
print(round(tab7[, c("week", "published", "pub_lo", "pub_hi", "pp_int", "pp_int_lo", "pp_int_hi", "rw1_int", "global_int")], 3))
cat("\nmean absolute difference from the published point estimates:\n")
print(round(sapply(names(fits), function(nm) mean(abs(tab7[[nm]] - tab7$published))), 3))
cat("The published estimates are the study authors' own probabilistic ML estimates, not a gold standard.\n")
matplot(weeks, tab7[, c("published", "pp_int", "rw1_int", "global_int")], type = "b", pch = 1:4, lty = 1,
        xlab = "week", ylab = "prevalence", main = "Weekly donor prevalence")
legend("topleft", c("published", "pp_int", "rw1_int", "global_int"), col = 1:4, pch = 1:4, lty = 1, bty = "n")

dev.off()
saveRDS(list(hold = hold_tab, ppc = ppc_tab, p4 = p4_tab, cmp = cmp, prev = tab7, obs_sd_p4 = obs_sd,
             rhat = sapply(fits, `[[`, "rhat")), file.path(out_dir, "workflow_results.rds"))
cat("\nFigures and results written to", out_dir, "\n")
