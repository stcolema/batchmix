# Semi-synthetic test on Stockholm ODs with KNOWN weekly prevalence.
# Donor items are real ODs resampled from the labelled controls: positives
# from the COVID-positive pool, negatives from the historical negatives of the
# same assay run, mixed at the published weekly prevalence (a realistic
# temporal shape with abrupt run-to-run changes). Weeks are small
# (n_per_week donors), so the weight prior matters. A few controls per
# class stay fixed as anchors, as in stockholm_serology.R.
#
#   Rscript stockholm_semisynthetic.R /path/to/seroprevalence-paper [reps] [n_per_week] [n_iter] [out.rds] [pos_quantile]
#
# pos_quantile (default 1): keep only the COVID-positive controls whose spike OD
# lies below this quantile of the pool. With 1 the classes are almost perfectly
# separable and every weight prior ties; 0.5 (the weaker half) makes donors
# resemble mild/waning infections, where the weight prior can matter.
suppressMessages({library(batchmix); library(data.table)})
args <- commandArgs(TRUE)
path <- if (length(args) >= 1) args[1] else "seroprevalence-paper"
reps <- if (length(args) >= 2) as.integer(args[2]) else 8L
n_w <- if (length(args) >= 3) as.integer(args[3]) else 30L
n_iter <- if (length(args) >= 4) as.integer(args[4]) else 2000L
out_file <- if (length(args) >= 5) args[5] else "stockholm_semisynthetic.rds"
pos_q <- if (length(args) >= 6) as.numeric(args[6]) else 1

e <- new.env(); load(file.path(path, "adjusted-data.RData"), e); m <- as.data.table(e$m)
est <- fread(file.path(path, "estimates.csv"))
m <- m[!duplicated(paste(group, type, Sample.ID))]
weeks <- c(14, 17:25, 30:34, 45:50)
run_of_week <- function(w) ifelse(w <= 21, "early", ifelse(w <= 25, "200702", ifelse(w <= 34, "200923", "201216")))
ref <- est[type == "Blood donors" & grepl("^Wk", Week)]; ref[, week := as.integer(sub("Wk", "", Week))]
p_true <- ref$estimate[match(weeks, ref$week)]
neg <- m[type == "Historical controls"]
neg[, run := ifelse(group %in% c("Blood Donors", "Pregnant volunteers"), "early", group)]
neg <- neg[run %in% run_of_week(weeks)] # run 201007 has no donor weeks
pos <- m[type == "COVID"]
pos <- pos[OD.spike <= quantile(OD.spike, pos_q)]
feat <- function(d) log10(cbind(pmax(d$OD.spike, 0.01), pmax(d$OD.rbd, 0.01)))
specs <- list(global = list("global"), pp = list("partial_pooling"),
              rw1 = list("gp", "rw1"), rw2 = list("gp", "rw2"), matern32 = list("gp", "matern32"))

one_rep <- function(seed) {
  set.seed(seed)
  pos_anchor_idx <- sample(nrow(pos), 40); pos_pool <- pos[-pos_anchor_idx]
  rows <- list(); truth <- numeric(length(weeks))
  for (i in seq_along(weeks)) {
    n_pos <- rbinom(1, n_w, p_true[i]); truth[i] <- n_pos / n_w
    r <- run_of_week(weeks[i])
    npool <- neg[run == r]
    rows[[i]] <- rbind(pos_pool[sample(.N, n_pos, TRUE)], npool[sample(.N, n_w - n_pos, TRUE)], fill = TRUE)
    rows[[i]][, `:=`(batch = i, anchor = 0L, lab = c(rep(1L, n_pos), rep(0L, n_w - n_pos)))]
  }
  don <- rbindlist(rows, fill = TRUE)
  # anchors: positives anywhere, negatives from the matching run
  a_pos <- pos[pos_anchor_idx]; a_pos[, `:=`(batch = sample(seq_along(weeks), .N, TRUE), anchor = 1L, lab = 1L)]
  a_neg <- neg[sample(.N, 80)]
  a_neg[, batch := vapply(run, function(r) sample(rep(which(run_of_week(weeks) == r), 2), 1), 1L)]
  a_neg[, `:=`(anchor = 1L, lab = 0L)]
  d <- rbindlist(list(don, a_pos, a_neg), fill = TRUE)
  X <- feat(d); fixed <- d$anchor
  init <- ifelse(fixed == 1, d$lab, as.integer(X[, 1] > log10(0.25)))
  sapply(names(specs), function(nm) {
    s <- specs[[nm]]
    a <- list(X, 3, n_iter, max(1, n_iter %/% 200), d$batch, "MVN", initial_labels = init, fixed = fixed,
              control = batchmixControl(n_burn = n_iter %/% 3), K_max = 2, batch_weight_prior = s[[1]])
    if (s[[1]] == "gp") a <- c(a, list(batch_coordinates = weeks, sample_gp_hyperparameters = TRUE, gp_kernel = s[[2]]))
    o <- suppressMessages(do.call(fitBatchMix, a))
    pb <- processMCMCChain(getBestChain(o), n_iter %/% 2)
    lab <- pb$pred - 1L
    est_p <- tapply(lab[d$anchor == 0], d$batch[d$anchor == 0], mean)
    c(mae = mean(abs(est_p - truth)), bias = mean(est_p - truth),
      rmse = sqrt(mean((est_p - truth)^2)),
      err_jump = mean(abs(est_p - truth)[c(11, 16)]))   # first week of runs 200923, 201216
  })
}
all <- lapply(seq_len(reps), function(i) { cat("rep", i, "\n"); one_rep(1000 + i) })
arr <- simplify2array(all)             # metric x model x rep
cat("\nMean (SE) over", reps, "replicates; truth = realised weekly fraction positive among donors\n")
for (met in dimnames(arr)[[1]]) {
  cat(met, ":", paste(sprintf("%s %.3f (%.3f)", dimnames(arr)[[2]], rowMeans(arr[met, , ]), apply(arr[met, , ], 1, sd) / sqrt(reps)), collapse = " | "), "\n")
}
cat("\nPaired MAE difference vs partial pooling (negative = better), mean (SE):\n")
for (nm in names(specs)[-2]) {
  dd <- arr["mae", nm, ] - arr["mae", "pp", ]
  cat(sprintf("%-9s %.4f (%.4f)\n", nm, mean(dd), sd(dd) / sqrt(reps)))
}
saveRDS(arr, out_file)
