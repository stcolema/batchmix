args <- commandArgs(TRUE)
is_t <- as.logical(args[1]); n_sc <- as.integer(args[2]); thin <- as.integer(args[3])
smS <- as.logical(args[4]); smM <- as.logical(args[5])
set.seed(as.integer(args[6])); Sys.setenv(PKG_CPPFLAGS = paste0("-I", normalizePath("../../../src")))
Rcpp::sourceCpp("geweke_mvn.cpp", verbose = FALSE, rebuild = FALSE)
t0 <- Sys.time()
o <- gewekeMVN(is_t, 400000, n_sc, thin, 0.7, 40, 0.5, 40, 40, smM, smS, as.integer(args[6]), as.logical(args[7]), as.logical(args[8]))
cat("elapsed", format(Sys.time() - t0), "\n")
cat("acc cov", round(o$cov_acc,3), " mu", round(o$mu_acc,3), " m", round(o$m_acc,3), " S", round(o$S_acc,3), if (is_t) c(" df", round(o$df_acc,3)), "\n")
names <- c("mu00","mu11","logdet0","cov0_00","corr0","m00","m11","S00","S11","N0","w0","mu00*m00","df0")
bm_se <- function(x, nb = 40) { b <- split(x, cut(seq_along(x), nb, labels = FALSE)); sd(vapply(b, mean, 0)) / sqrt(nb) }
cat(sprintf("%-9s %10s %10s %8s | quantile indicators z (5,25,50,75,95%%)\n", "stat", "MC mean", "SC mean", "z(mean)"))
for (j in seq_len(ncol(o$mc))) {
  a <- o$mc[, j]; b <- o$sc[, j]
  se_mc <- sd(a) / sqrt(length(a)); se_sc <- bm_se(b)
  zmean <- (mean(b) - mean(a)) / sqrt(se_mc^2 + se_sc^2)
  zq <- sapply(c(.05, .25, .5, .75, .95), function(p) { q <- quantile(a, p); ind <- as.numeric(b <= q)
    (mean(ind) - p) / max(bm_se(ind), sqrt(p*(1-p)/length(ind))) })
  cat(sprintf("%-9s %10.4f %10.4f %8.2f | %s\n", names[j], mean(a), mean(b), zmean, paste(sprintf("%5.1f", zq), collapse = " ")))
}
