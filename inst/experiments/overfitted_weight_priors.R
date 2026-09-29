# Over-fitted mixtures under different priors on the mixture weights.
#
# Data: 3 well-separated 2-d clusters in 3 batches (N = 450), fitted with
# K_max = 8, i.e. 5 more components than exist. For each weight prior this
# reports the number of occupied components after burn-in (kocc_*) and the
# fraction of items outside the 3 largest clusters (extra_frac); a prior that
# empties surplus components gives kocc near 3 and extra_frac near 0. Four
# seeds, one chain of 4000 iterations each, so treat it as a demonstration of
# the direction and rough size of the effect, not a precise estimate.
#
# Result recorded when this was written (means over the four seeds; partial
# pooling with K exchangeable logits and no reference class):
#   global, alpha = 1/K_max   kocc 4.08  extra_frac 0.063
#   global, alpha = 1         kocc 7.27  extra_frac 0.226
#   partial pooling, sd 10    kocc 3.81  extra_frac 0.036   (the default)
#   partial pooling, sd 2.5   kocc 4.83  extra_frac 0.048
#   partial pooling, sd 1     kocc 7.85  extra_frac 0.269
# (An earlier parameterisation against a reference class gave 4.21, 6.46 and
# 7.90 for sd 10, 2.5 and 1: it was not exchangeable over clusters and less
# sparse at moderate sd.)
# The diffuse default logit prior behaves like a sparse Dirichlet; a much
# sharper one does not. Rousseau & Mengersen (2011, JRSS-B 73(5))
# prove the emptying behaviour for Dirichlet priors with concentration below
# half the component parameter dimension; nothing here is a theorem for the
# logistic-normal prior.

library(batchmix)
library(parallel)
K_true <- 3; K_max <- 8; N <- 450; B <- 3; P <- 2
sim <- function(seed) {
  set.seed(seed)
  s <- generateBatchData(N, P, c(0, 5, 10), rep(1, K_true), rnorm(B, 0.2, 0.15), rep(1.2, B),
                         rep(1/K_true, K_true), rep(1/B, B), type = "MVN")
  list(X = s$observed_data, b = s$batch_IDs)
}
configs <- list(
  global_alpha_1overK = list(batch_weight_prior = "global", alpha = 1 / K_max),
  global_alpha_1      = list(batch_weight_prior = "global", alpha = 1),
  pp_sd10             = list(batch_weight_prior = "partial_pooling", pp_mu_prior_sd = 10),
  pp_sd2.5            = list(batch_weight_prior = "partial_pooling", pp_mu_prior_sd = 2.5),
  pp_sd1              = list(batch_weight_prior = "partial_pooling", pp_mu_prior_sd = 1)
)
jobs <- expand.grid(cfg = names(configs), seed = 1:4, stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(jobs)), function(i) { tryCatch({
  d <- sim(100 + jobs$seed[i]); cf <- configs[[jobs$cfg[i]]]
  set.seed(jobs$seed[i])
  out <- suppressWarnings(do.call(runBatchMix, c(list(d$X, n_iter = 4000, thin = 20, batch_vec = d$b, type = "MVN",
                    K_max = K_max, verbose = FALSE), cf)))
  lab <- out$samples[(nrow(out$samples)/2 + 1):nrow(out$samples), , drop = FALSE]
  kocc <- apply(lab, 1, function(l) length(unique(l)))
  # size of clusters beyond the 3 largest, as a fraction of N
  extra <- apply(lab, 1, function(l) { t <- sort(table(l), decreasing = TRUE); if (length(t) > 3) sum(t[-(1:3)]) / N else 0 })
  c(kocc_mean = mean(kocc), kocc_min = min(kocc), extra_frac = mean(extra))
}, error = function(e) c(kocc_mean = NA, kocc_min = NA, extra_frac = NA, msg = conditionMessage(e)))}, mc.cores = 4)
tab <- cbind(jobs, do.call(rbind, res))
print(aggregate(cbind(kocc_mean, kocc_min, extra_frac) ~ cfg, tab, mean), digits = 3)
print(tab, digits = 3)
