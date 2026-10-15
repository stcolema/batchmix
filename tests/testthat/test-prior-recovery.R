# Prior recovery for every update block (cf. Talts et al., 2018, arXiv:1804.06788,
# on checking a sampler against its own prior). The chains in
# src/priorOnlyChains.cpp swap the data likelihood for zero, so each block's
# Metropolis or Gibbs kernel is left with its prior terms, Jacobians and proposal
# corrections only, and the chain's stationary law is the prior, known in closed
# form. A kernel that is wrong in any of these shifts a marginal away from it.
#
# Checks are on indicators of {x <= q}, with q the exact quantiles of the prior,
# not on moments: several priors here (inverse-Wishart, inverse-gamma) have
# infinite variance, which would make a standard error from the chain
# meaningless. The standard error is that of batch means, floored at the
# independent-draws value.

bm_se <- function(x, n_batches = 40) {
  batches <- split(x, cut(seq_along(x), n_batches, labels = FALSE))
  stats::sd(vapply(batches, mean, numeric(1))) / sqrt(n_batches)
}

expect_follows <- function(x, qfun, label, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), z_max = 4.5) {
  for (p in probs) {
    ind <- as.numeric(x <= qfun(p))
    se <- max(bm_se(ind), sqrt(p * (1 - p) / length(ind)))
    expect_lt(abs(mean(ind) - p) / se, z_max, label = sprintf("%s, P(X <= q%g)", label, p))
  }
}

q_invgamma <- function(p, shape, scale) 1 / stats::qgamma(1 - p, shape = shape, rate = scale)

# Quantiles of S - S_loc when the concentration rho is sampled: with rho - 2 ~ Gamma(a_s, b_s)
# and theta = s_scale_prior_mean (rho - 1), a mixture of InvGamma(rho, theta) over rho
s_scale_mixture_quantile <- function(out) {
  cdf <- function(x) {
    stats::integrate(function(r) {
      vapply(r, function(ri) {
        stats::pgamma(1 / x, shape = ri, rate = out$s_scale_prior_mean * (ri - 1), lower.tail = FALSE) *
          stats::dgamma(ri - 2, out$a_s, out$b_s)
      }, numeric(1))
    }, lower = 2, upper = Inf, rel.tol = 1e-8)$value
  }
  function(p) stats::uniroot(function(x) cdf(x) - p, c(1e-6, 1e4), tol = 1e-10)$root
}

# Occupied clusters 0 and 1 (every other cluster empty), every batch occupied
prior_run <- function(type = "MVN", K = 2, B = 2, N = 8, P = 2, n_iter = 100000, seed = 1,
                      sample_m_scale = FALSE, sample_s_scale = FALSE, rho = 3, theta = 1,
                      m_scale = 0.1, update_allocation = FALSE, concentration = rep(1, K)) {
  set.seed(seed)
  X <- matrix(stats::rnorm(N * P), N, P)
  priorOnlyChain(
    type = type, X = X, K = K, B = B,
    labels = rep_len(0:1, N), batch_vec = rep_len(0:(B - 1), N), n_iter = n_iter,
    mu_pw = 0.5, cov_pw = 50, m_pw = 0.5, S_pw = 50, t_df_pw = 30, m_scale = m_scale,
    rho = rho, theta = theta, sample_m_scale = sample_m_scale, sample_s_scale = sample_s_scale,
    s_scale_pw = 0.5, concentration = concentration, update_allocation = update_allocation
  )
}

# Checks the cluster-level blocks of cluster `k` (0-based): its covariance and mean
check_cluster_prior <- function(out, k = 0, label = "") {
  P <- nrow(out$xi)
  Psi <- out$scale
  S11 <- out$cov[1, k * P + 1, ]
  # Sigma ~ InvWishart(Psi, nu): the (1, 1) element is InvGamma((nu - P + 1) / 2, Psi_11 / 2)
  expect_follows(S11, function(p) q_invgamma(p, (out$nu - P + 1) / 2, Psi[1, 1] / 2),
                 paste(label, "covariance (1, 1) element"))
  # The log determinant: E log|Sigma| = log|Psi| - P log 2 - sum_i digamma((nu + 1 - i) / 2)
  # (a reference value, with the batch-means standard error of the mean)
  logdet <- vapply(seq_len(dim(out$cov)[3]), function(i) {
    as.numeric(determinant(out$cov[, k * P + seq_len(P), i], logarithm = TRUE)$modulus)
  }, numeric(1))
  expected <- as.numeric(determinant(Psi, logarithm = TRUE)$modulus) - P * log(2) -
    sum(digamma((out$nu + 1 - seq_len(P)) / 2))
  expect_lt(abs(mean(logdet) - expected) / bm_se(logdet), 4.5, label = paste(label, "E log|Sigma|"))
  # kappa (mu - xi)' Sigma^-1 (mu - xi) ~ chi^2_P, whatever Sigma
  stat <- vapply(seq_len(dim(out$cov)[3]), function(i) {
    d <- out$mu[, k + 1, i] - out$xi
    out$kappa * drop(t(d) %*% solve(out$cov[, k * P + seq_len(P), i], d))
  }, numeric(1))
  expect_follows(stat, function(p) stats::qchisq(p, P), paste(label, "mean given covariance"))
}

test_that("MVN: cluster mean and covariance of an occupied cluster follow the prior", {
  out <- prior_run("MVN", seed = 101)
  expect_gt(out$cov_count[1], 1000)
  expect_gt(out$mu_count[1], 1000)
  check_cluster_prior(out, 0, "MVN occupied")
})

test_that("MVN: batch shift follows N(0, delta_2 lambda_2), and lambda_2 its Inverse-Gamma hyperprior", {
  out <- prior_run("MVN", seed = 102, m_scale = 0.1)
  sd_m <- sqrt(out$delta_2 * 0.1)
  expect_gt(out$m_count[1], 1000)
  expect_follows(out$m[1, 1, ] / sd_m, stats::qnorm, "batch shift, fixed lambda_2")
  expect_follows(out$m[2, 2, ] / sd_m, stats::qnorm, "batch shift, second cell")

  out <- prior_run("MVN", seed = 103, sample_m_scale = TRUE)
  # lambda_2 ~ InvGamma(a, b); m | lambda_2 ~ N(0, delta_2 lambda_2), hence m is a scaled t
  # with 2a degrees of freedom and scale sqrt(delta_2 b / a)
  expect_follows(out$lambda_2, function(p) q_invgamma(p, out$a, out$b), "lambda_2")
  expect_follows(out$m[1, 1, ] / sqrt(out$delta_2 * out$b / out$a),
                 function(p) stats::qt(p, df = 2 * out$a), "batch shift, sampled lambda_2")
})

test_that("MVN: batch scale follows S_loc + InvGamma(rho, theta), also when its concentration is sampled", {
  out <- prior_run("MVN", seed = 104, rho = 3, theta = 1)
  expect_gt(out$S_count[1], 1000)
  expect_follows(out$S[1, 1, ] - out$S_loc, function(p) q_invgamma(p, 3, 1), "batch scale, fixed rho")
  expect_follows(out$S[2, 2, ] - out$S_loc, function(p) q_invgamma(p, 3, 1), "batch scale, second cell")
  expect_true(all(out$rho == 3))

  # rho - 2 ~ Gamma(a_s, b_s), and S - S_loc a mixture over it (see s_scale_mixture_quantile())
  out <- prior_run("MVN", seed = 105, rho = 3, theta = 1, sample_s_scale = TRUE)
  expect_gt(out$s_scale_count, 1000)
  expect_follows(out$rho - 2, function(p) stats::qgamma(p, out$a_s, out$b_s), "rho")
  expect_follows(out$S[1, 1, ] - out$S_loc, s_scale_mixture_quantile(out), "batch scale, sampled rho")
})

test_that("MVT: cluster, batch and degrees-of-freedom blocks follow the prior", {
  # K = 3: clusters 0 and 1 are occupied (Metropolis-Hastings), cluster 2 is empty (forced prior draws)
  out <- prior_run("MVT", K = 3, seed = 106, sample_m_scale = TRUE, sample_s_scale = TRUE)
  expect_gt(out$t_df_count[1], 1000)
  check_cluster_prior(out, 0, "MVT occupied")
  check_cluster_prior(out, 2, "MVT empty")
  # df - t_loc ~ Gamma(psi, chi)
  expect_follows(out$t_df[, 1] - out$t_loc, function(p) stats::qgamma(p, out$psi, out$chi), "t df, occupied")
  expect_follows(out$t_df[, 3] - out$t_loc, function(p) stats::qgamma(p, out$psi, out$chi), "t df, empty")
  expect_follows(out$S[1, 1, ] - out$S_loc, s_scale_mixture_quantile(out), "MVT batch scale, sampled rho")
  expect_follows(out$lambda_2, function(p) q_invgamma(p, out$a, out$b), "MVT lambda_2")
})

test_that("empty clusters take exact prior draws while the labels move", {
  # K = 3 clusters, so that a cluster is sometimes empty, sometimes occupied
  out <- prior_run("MVN", K = 3, B = 2, N = 4, seed = 107, n_iter = 100000, update_allocation = TRUE)
  keep <- 2001:100000
  one <- out
  one$cov <- out$cov[, , keep]
  one$mu <- out$mu[, , keep]
  check_cluster_prior(one, 0, "MVN, labels moving: cluster 0")
  check_cluster_prior(one, 2, "MVN, labels moving: cluster 2")
})

test_that("labels and weights follow the Dirichlet-multinomial prior", {
  # N = 5 items, K = 3 clusters, weights ~ Dirichlet(alpha), labels | weights iid categorical
  alpha <- c(0.7, 1.4, 2.1)
  K <- 3
  N <- 5
  out <- prior_run("MVN", K = K, B = 1, N = N, P = 1, seed = 108, n_iter = 200000,
                   update_allocation = TRUE, concentration = alpha)
  keep <- 2001:200000
  counts <- t(apply(out$labels[keep, , drop = FALSE] + 1, 1, tabulate, nbins = K)) # labels are 0-based
  A <- sum(alpha)
  # P(n_1, ..., n_K) = N! / prod(n_k!) * Gamma(A) / Gamma(A + N) * prod Gamma(alpha_k + n_k) / Gamma(alpha_k)
  compositions <- as.matrix(expand.grid(0:N, 0:N, 0:N))
  compositions <- compositions[rowSums(compositions) == N, , drop = FALSE]
  log_p <- apply(compositions, 1, function(n) {
    lgamma(N + 1) - sum(lgamma(n + 1)) + lgamma(A) - lgamma(A + N) + sum(lgamma(alpha + n) - lgamma(alpha))
  })
  expect_equal(sum(exp(log_p)), 1, tolerance = 1e-12)
  worst <- 0
  for (j in seq_len(nrow(compositions))) {
    ind <- as.numeric(counts[, 1] == compositions[j, 1] & counts[, 2] == compositions[j, 2] &
                        counts[, 3] == compositions[j, 3])
    p <- exp(log_p[j])
    se <- max(bm_se(ind), sqrt(p * (1 - p) / length(ind)))
    worst <- max(worst, abs(mean(ind) - p) / se)
  }
  expect_lt(worst, 4.5, label = "largest z over the 21 compositions")
  # The weights: w_k ~ Beta(alpha_k, A - alpha_k)
  for (k in 1:K) {
    expect_follows(out$w[keep, k], function(p) stats::qbeta(p, alpha[k], A - alpha[k]), paste("weight", k))
  }
})

interaction_checks <- function(type, K, B, seed) {
  set.seed(seed)
  a <- 3
  b <- 2
  out <- priorOnlyInteractionChain(type, P = 2, K = K, B = B, n_iter = 100000, gamma_pw = 0.5,
                                   a_gamma = a, b_gamma = b)
  keep <- 2001:100000
  expect_gt(out$gamma_count[1], 1000)
  # tau2 ~ InvGamma(a_gamma, b_gamma)
  expect_follows(out$tau2[keep, 1], function(p) q_invgamma(p, a, b), paste(type, "tau2_interaction"))
  # Given tau2, gamma[p, , ] is Gaussian on the (K - 1)(B - 1)-dimensional sum-to-zero
  # subspace, so its sum of squares over tau2 is chi^2 with that many degrees of freedom...
  d <- (K - 1) * (B - 1)
  expect_follows(out$gamma_ss[keep, 1] / out$tau2[keep, 1], function(p) stats::qchisq(p, d),
                 paste(type, "sum of squares of gamma over tau2"))
  # ... and each entry has variance tau2 (1 - 1 / K)(1 - 1 / B)
  z <- out$gamma[1, 1, keep] / sqrt(out$tau2[keep, 1] * (1 - 1 / K) * (1 - 1 / B))
  expect_follows(z, stats::qnorm, paste(type, "entry of gamma"))
  # every row and column sums to zero
  expect_lt(max(abs(rowSums(out$gamma[, , 5000]))), 1e-10)
  expect_lt(max(abs(colSums(out$gamma[, , 5000]))), 1e-10)
}

test_that("the interaction term and its variance follow the prior", {
  interaction_checks("MVN", K = 3, B = 4, seed = 109)
  interaction_checks("MVT", K = 2, B = 3, seed = 110)
})

test_that("the Gaussian process hyperparameters follow their priors", {
  coords <- c(0, 0.15, 0.4, 0.7, 1)
  for (kernel_type in 0:2) {
    set.seed(120 + kernel_type)
    out <- gpWeightChain(K = 2, B = 5, n_iter = 100000, coords = coords, kernel_type = kernel_type,
                         tau2 = 1, length_scale = 1, level_sd = 3, eta_pw = 0.5,
                         sample_hyperparameters = TRUE, hyper_pw = 0.5)
    expect_gt(out$gp_hyperparameter_count, 1000)
    keep <- 2001:100000
    expect_follows(out$gp_tau2[keep], function(p) q_invgamma(p, out$gp_tau2_prior_shape, out$gp_tau2_prior_rate),
                   paste("GP tau2, kernel", kernel_type))
    if (kernel_type == 0) {
      expect_follows(out$gp_length_scale[keep],
                     function(p) q_invgamma(p, out$gp_length_scale_prior_shape, out$gp_length_scale_prior_rate),
                     "GP length scale, Matern-3/2")
    }
  }
})
