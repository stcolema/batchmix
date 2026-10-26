args <- commandArgs(TRUE)
mixed <- as.logical(args[1]); n_sc <- as.integer(args[2]); thin <- as.integer(args[3]); set.seed(as.integer(args[4]))
Sys.setenv(PKG_CPPFLAGS = paste0("-I", normalizePath("../../../src")))
Rcpp::sourceCpp("geweke_sep.cpp", verbose = FALSE)
t0 <- Sys.time()
o <- gewekeSep(mixed, 400000, n_sc, thin, 0.4, 40, 0.7, 0.5, 40, 2)
cat("elapsed", format(Sys.time() - t0), "\n")
cat("acc r", round(o$r_acc,3), " sigma", round(o$sigma_acc,3), " mu", round(o$mu_acc,3), " m", round(o$m_acc,3), " S", round(o$S_acc,3), "\n")
names <- c("mu00","mu11","logdet0","cov0_00","R01","m00","m11","S00","S11","w0","sigma_1(k0)")
bm_se <- function(x, nb = 40) { b <- split(x, cut(seq_along(x), nb, labels = FALSE)); sd(vapply(b, mean, 0)) / sqrt(nb) }
cat(sprintf("%-11s %9s %9s %8s | quantile z (5,25,50,75,95%%)\n", "stat", "MC mean", "SC mean", "z(mean)"))
for (j in seq_len(ncol(o$mc))) {
  a <- o$mc[, j]; b <- o$sc[, j]
  se_mc <- sd(a) / sqrt(length(a)); se_sc <- bm_se(b)
  zmean <- (mean(b) - mean(a)) / sqrt(se_mc^2 + se_sc^2)
  zq <- sapply(c(.05, .25, .5, .75, .95), function(p) { q <- quantile(a, p); ind <- as.numeric(b <= q)
    (mean(ind) - p) / max(bm_se(ind), sqrt(p*(1-p)/length(ind))) })
  cat(sprintf("%-11s %9.4f %9.4f %8.2f | %s\n", names[j], mean(a), mean(b), zmean, paste(sprintf("%5.1f", zq), collapse = " ")))
}
