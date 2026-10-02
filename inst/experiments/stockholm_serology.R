# Stockholm SARS-CoV-2 serology (De Castro et al.; chr1swallace/seroprevalence-paper):
# blood donors, 100 per week over 21 weeks, raw ODs for spike and RBD, plus
# historical negative and COVID-positive controls. Compares batch-weight priors.
#
#   Rscript stockholm_serology.R /path/to/seroprevalence-paper [n_iter] [out.rds]
#
# Design (everything below is a choice, not a given):
#  * batch = week (21 batches; weeks nest inside four assay runs, so week
#    and run effects are confounded - the plate/run is not in the file).
#  * features: log10 of raw OD (spike, RBD), not the plate-adjusted columns.
#  * controls have no week: historical negatives are assigned uniformly at
#    random to weeks of the run they were measured in (run 201007 has no
#    donor weeks and is dropped), COVID-positives to any week. Half of the
#    controls keep their label (fixed); the other half are unlabelled and
#    used to score classification. Controls count in each batch's class
#    counts, which biases the weight prior toward the controls' mix, so the
#    prevalence reported is the fraction of DONOR samples allocated to the
#    positive class, not the batch weight.
suppressMessages({library(batchmix); library(data.table)})
args <- commandArgs(TRUE)
path <- if (length(args) >= 1) args[1] else "seroprevalence-paper"
n_iter <- if (length(args) >= 2) as.integer(args[2]) else 3000L
out_file <- if (length(args) >= 3) args[3] else "stockholm_results.rds"
set.seed(2024)

e <- new.env(); load(file.path(path, "adjusted-data.RData"), e); m <- as.data.table(e$m)
est <- fread(file.path(path, "estimates.csv"))
m[, week := suppressWarnings(as.integer(ifelse(grepl("Wk", Sample.ID), sub(".*Wk([0-9]+).*", "\\1", Sample.ID), NA)))]
m <- m[!duplicated(paste(group, type, Sample.ID))]

donors <- m[type == "Blood donors"]
run_of_week <- function(w) ifelse(w <= 21, "early", ifelse(w <= 25, "200702", ifelse(w <= 34, "200923", "201216")))
donors[, run := run_of_week(week)]
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
X <- log10(as.matrix(d[, .(pmax(OD.spike, 0.01), pmax(OD.rbd, 0.01))]))
batch <- match(d$week, weeks)
fixed <- as.integer(d$role == "fixed")
init <- ifelse(fixed == 1, d$truth, as.integer(X[, 1] > log10(0.25)))
cat("N =", nrow(d), "| donors", sum(d$role == "donor"), "| fixed", sum(fixed), "| holdout", sum(d$role == "holdout"), "\n")

fit <- function(prior, kernel = NULL) {
  a <- list(X, 4, n_iter, max(1, n_iter %/% 200), batch, "MVN", initial_labels = init, fixed = fixed,
            control = batchmixControl(n_burn = n_iter %/% 3), K_max = 2, batch_weight_prior = prior)
  if (prior == "gp") a <- c(a, list(batch_coordinates = weeks, sample_gp_hyperparameters = TRUE, gp_kernel = kernel))
  o <- do.call(fitBatchMix, a)
  p <- processMCMCChains(o, n_iter %/% 2)
  best <- getBestChain(o); pb <- processMCMCChain(best, n_iter %/% 2)
  lab <- pb$pred - 1L
  is_d <- d$role == "donor"; is_h <- d$role == "holdout"
  # item-level posterior probability of the positive class (donors): average of allocation draws
  alloc <- if (!is.null(pb$allocation)) pb$allocation else NULL
  prev <- tapply(lab[is_d], d$week[is_d], mean)
  list(prev = prev, sens = mean(lab[is_h & d$truth == 1] == 1), spec = mean(lab[is_h & d$truth == 0] == 0),
       rhat = attr(o, "convergence")$rhat, w = if (!is.null(pb$w_batch_est)) pb$w_batch_est[, 2] else NULL)
}
specs <- list(global = list("global"), pp = list("partial_pooling"),
              rw1 = list("gp", "rw1"), rw2 = list("gp", "rw2"), matern32 = list("gp", "matern32"))
res <- lapply(specs, function(s) { t0 <- Sys.time(); r <- do.call(fit, s); r$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs")); r })
ref <- est[type == "Blood donors" & grepl("^Wk", Week)]; ref[, week := as.integer(sub("Wk", "", Week))]
ref <- ref[match(weeks, week)]
tab <- do.call(rbind, lapply(names(res), function(n) {
  r <- res[[n]]; pr <- r$prev
  data.frame(model = n, sens = r$sens, spec = r$spec,
             MAE_vs_published = mean(abs(pr - ref$estimate), na.rm = TRUE),
             in_published_CI = mean(pr >= ref$lower.ci & pr <= ref$upper.ci, na.rm = TRUE),
             rhat = r$rhat, secs = r$secs)
}))
print(tab, digits = 3)
print(round(cbind(week = weeks, published = ref$estimate, sapply(res, function(r) as.numeric(r$prev))), 3))
saveRDS(list(tab = tab, res = res, weeks = weeks, ref = ref), out_file)
