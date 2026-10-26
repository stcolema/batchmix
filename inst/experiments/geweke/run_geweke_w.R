args <- commandArgs(TRUE)
type <- as.integer(args[1]); kernel <- as.integer(args[2]); hyper <- as.logical(args[3]); n_sc <- as.integer(args[4]); set.seed(as.integer(args[5]))
Sys.setenv(PKG_CPPFLAGS = paste0("-I", normalizePath("../../../src")))
Rcpp::sourceCpp("geweke_w.cpp", verbose = FALSE)
t0 <- Sys.time()
o <- gewekeW(type, kernel, hyper, 400000, n_sc, 5, 0.6, 0.5, if (type == 1) 2 else 3, 8)
cat("type", type, "kernel", kernel, "hyper", hyper, " elapsed", format(Sys.time() - t0), " eta acc", round(o$eta_acc, 3), " hyper acc", round(o$hyper_acc, 3), "\n")
names <- c("d(0:1)","d(1:0,2)","d(B:1-0)","mu0/beta0","log tau2_0","log tau2_1/ell","w[0,0]","w[B,2]")
bm_se <- function(x, nb = 40) { b <- split(x, cut(seq_along(x), nb, labels = FALSE)); sd(vapply(b, mean, 0)) / sqrt(nb) }
for (j in seq_len(ncol(o$mc))) {
  a <- o$mc[, j]; b <- o$sc[, j]
  se_mc <- sd(a) / sqrt(length(a)); se_sc <- bm_se(b)
  zmean <- (mean(b) - mean(a)) / sqrt(se_mc^2 + se_sc^2)
  zv <- (log(var(b)) - log(var(a))) / sqrt(2 * (bm_se((b-mean(b))^2) / var(b))^2 + 2 * (se_mc*sqrt(length(a))/sd(a))^2 * 0 + 1e-12)
  zq <- sapply(c(.05, .25, .5, .75, .95), function(p) { q <- quantile(a, p); ind <- as.numeric(b <= q)
    (mean(ind) - p) / max(bm_se(ind), sqrt(p*(1-p)/length(ind))) })
  cat(sprintf("  %-13s MC %8.4f SC %8.4f z(mean) %5.2f | var MC %7.4f SC %7.4f | q-z %s\n", names[j], mean(a), mean(b), zmean, var(a), var(b), paste(sprintf("%5.1f", zq), collapse = " ")))
}
