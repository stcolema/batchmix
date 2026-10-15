# The label posterior of a tiny model against exact enumeration.
#
# One feature (P = 1), one batch (B = 1), K = 2 clusters, global weights, fixed
# batch-shift variance (lambda_2) and fixed scale prior (rho, theta). The model is
#   w ~ Dirichlet(alpha, alpha),   z_n | w ~ Categorical(w),
#   sigma2_k ~ InvGamma(nu / 2, psi / 2),   mu_k | sigma2_k ~ N(xi, sigma2_k / kappa),
#   m ~ N(0, delta_2 lambda_2),   S - 1 ~ InvGamma(rho, theta),
#   x_n | z_n = k ~ N(mu_k + m, S sigma2_k),
# with xi, psi = var(x) / K^2 and delta_2 = var(x) set from the data as in the sampler.
#
# Given (m, S), each cluster's (mu, sigma2) integrates out in closed form (a Student-t
# marginal), and the weights integrate out to a Dirichlet-multinomial. What remains,
# an integral over (m, S), is done by Gauss quadrature, so p(z | x) is known for each
# of the 2^N labellings of N items. The uncollapsed sampler, which updates the
# parameters and the labels in turn, must reproduce it: its likelihood, every
# Metropolis kernel and the label and weight updates are all involved.

bm_se_exact <- function(x, n_batches = 40) {
  batches <- split(x, cut(seq_along(x), n_batches, labels = FALSE))
  stats::sd(vapply(batches, mean, numeric(1))) / sqrt(n_batches)
}

# Gauss-Legendre and Gauss-Hermite rules by the Golub-Welsch algorithm
gauss_rule <- function(n, type = c("legendre", "hermite")) {
  type <- match.arg(type)
  i <- seq_len(n - 1)
  off <- if (type == "legendre") i / sqrt(4 * i^2 - 1) else sqrt(i)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- off
  J[cbind(i + 1, i)] <- off
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  # nodes, and weights that sum to one (the expectation under U(-1, 1) or N(0, 1))
  list(x = e$values[ord], w = e$vectors[1, ord]^2)
}

# log of the marginal likelihood of the items `x` in one cluster, given the batch shift
# m and the batch scale S (vectors over a grid), with (mu, sigma2) integrated out
log_cluster_marginal <- function(x, m, S, xi, kappa, a0, b0) {
  n <- length(x)
  if (n == 0) return(matrix(0, length(m), length(S)))
  d <- outer(x, m, function(xx, mm) xx - mm) - xi # n x length(m)
  sd <- colSums(d)
  sd2 <- colSums(d^2)
  # q = d' Sigma0^-1 d with Sigma0 = S I + (1 / kappa) 11', by Sherman-Morrison
  q <- (matrix(sd2, length(m), length(S)) -
          outer(sd^2, 1 / kappa / (S + n / kappa), `*`)) / rep(S, each = length(m))
  logS <- matrix(log(S), length(m), length(S), byrow = TRUE)
  -0.5 * n * log(2 * pi) - 0.5 * ((n - 1) * logS + log(matrix(S + n / kappa, length(m), length(S), byrow = TRUE))) +
    a0 * log(b0) - lgamma(a0) + lgamma(a0 + n / 2) - (a0 + n / 2) * log(b0 + q / 2)
}

exact_label_posterior <- function(x, alpha = 1, m_scale = 0.1, rho = 3, theta = 1, kappa = 0.01,
                                  n_m = 120, n_s = 400) {
  K <- 2
  N <- length(x)
  v <- stats::var(x)
  xi <- mean(x)
  nu <- 3 # P + 2
  a0 <- nu / 2
  b0 <- (v / K^2) / 2
  # m ~ N(0, delta_2 lambda_2): Gauss-Hermite; S - 1 ~ InvGamma(rho, theta): Gauss-Legendre in its quantile
  gh <- gauss_rule(n_m, "hermite")
  m <- sqrt(v * m_scale) * gh$x
  gl <- gauss_rule(n_s, "legendre")
  u <- (gl$x + 1) / 2
  S <- 1 + 1 / stats::qgamma(1 - u, shape = rho, rate = theta)
  log_w <- outer(log(gh$w), log(gl$w), `+`)
  grid <- as.matrix(expand.grid(rep(list(0:1), N)))
  A <- K * alpha
  log_post <- apply(grid, 1, function(z) {
    n_k <- c(sum(z == 0), sum(z == 1))
    log_pz <- lgamma(A) - lgamma(A + N) + sum(lgamma(alpha + n_k) - lgamma(alpha))
    ll <- log_cluster_marginal(x[z == 0], m, S, xi, kappa, a0, b0) +
      log_cluster_marginal(x[z == 1], m, S, xi, kappa, a0, b0) + log_w
    mx <- max(ll)
    log_pz + mx + log(sum(exp(ll - mx)))
  })
  p <- exp(log_post - max(log_post))
  structure(p / sum(p), configs = grid)
}

test_that("the closed-form cluster marginal matches numerical integration", {
  x <- c(0.4, -0.3, 1.1)
  m <- 0.2
  S <- 1.7
  xi <- 0.1
  kappa <- 0.5
  a0 <- 1.5
  b0 <- 0.8
  closed <- log_cluster_marginal(x, m, S, xi, kappa, a0, b0)
  # integrate over sigma2 (InvGamma) and mu (Gaussian), numerically
  inner <- function(sigma2) {
    vapply(sigma2, function(s2) {
      stats::integrate(function(mu) {
        vapply(mu, function(mm) {
          prod(stats::dnorm(x, mm + m, sqrt(S * s2))) * stats::dnorm(mm, xi, sqrt(s2 / kappa))
        }, numeric(1))
      }, -Inf, Inf, rel.tol = 1e-10)$value * exp(a0 * log(b0) - lgamma(a0) - (a0 + 1) * log(s2) - b0 / s2)
    }, numeric(1))
  }
  numerical <- log(stats::integrate(inner, 0, Inf, rel.tol = 1e-8, subdivisions = 500)$value)
  expect_equal(drop(closed), numerical, tolerance = 1e-6)
  # an empty cluster carries no information
  expect_equal(drop(log_cluster_marginal(numeric(0), m, S, xi, kappa, a0, b0)), 0)
})

test_that("the quadrature is accurate: refining the grids does not change the posterior", {
  x <- c(-1.3, -0.6, 0.2, 0.7, 1.6)
  # the grid used below against one more than twice as fine: the error is about 4e-7, orders
  # below the Monte Carlo resolution (about 1e-3) of the comparison with the sampler
  used <- exact_label_posterior(x)
  fine <- exact_label_posterior(x, n_m = 320, n_s = 1200)
  expect_lt(max(abs(used - fine)), 2e-6)
  expect_equal(sum(fine), 1)
  # labels are exchangeable: swapping the two cluster names leaves the posterior unchanged
  configs <- attr(fine, "configs")
  swapped <- match(apply(1 - configs, 1, paste, collapse = ""), apply(configs, 1, paste, collapse = ""))
  expect_equal(as.numeric(fine), as.numeric(fine)[swapped], tolerance = 1e-8)
})

test_that("the sampler's label posterior matches exact enumeration", {
  x <- c(-1.3, -0.6, 0.2, 0.7, 1.6)
  exact <- exact_label_posterior(x, alpha = 1, m_scale = 0.1, rho = 3, theta = 1)
  configs <- attr(exact, "configs")
  # not concentrated: several labellings carry real mass, so the comparison has something to resolve
  expect_gt(sum(exact > 0.01), 6)

  set.seed(2024)
  fit <- runBatchMix(matrix(x, ncol = 1),
    n_iter = 600000, thin = 3, batch_vec = rep(1, length(x)), type = "MVN", K_max = 2, alpha = 1,
    m_scale = 0.1, rho = 3, theta = 1, sample_s_scale = FALSE, batch_weight_prior = "global",
    control = batchmixControl(auto_tune = FALSE), verbose = FALSE
  )
  labels <- fit$samples[-(1:2000), , drop = FALSE] # 0-based
  key <- apply(labels, 1, paste, collapse = "")
  worst <- 0
  for (j in seq_len(nrow(configs))) {
    ind <- as.numeric(key == paste(configs[j, ], collapse = ""))
    p <- exact[j]
    se <- max(bm_se_exact(ind), sqrt(p * (1 - p) / length(ind)))
    worst <- max(worst, abs(mean(ind) - p) / se)
  }
  expect_lt(worst, 4.5, label = "largest z over the 32 labellings")
  # a summary with its own resolution: the posterior mean number of items in cluster 0
  # (the label names are exchangeable, so use the size of the smaller cluster)
  smaller <- pmin(rowSums(labels), ncol(labels) - rowSums(labels))
  exact_smaller <- pmin(rowSums(configs), ncol(configs) - rowSums(configs))
  expect_lt(abs(mean(smaller) - sum(exact * exact_smaller)) / bm_se_exact(smaller), 4.5,
            label = "mean size of the smaller cluster")
})
