#
# Unit tests for the log-likelihood functions implemented in C++.
#
context("Likelihood comparison")
library(batchmix)

test_that("Gamma log-likelihood ratios.", {
  x1 <- 3
  x2 <- 4
  shape <- 3
  rate <- 1

  ratio1 <- gammaLogLikelihood(x1, shape, rate) - gammaLogLikelihood(x2, shape, rate)
  ratio2 <- dgamma(x1, shape, rate, log = T) - dgamma(x2, shape, rate, log = T)

  expect_equal(ratio1, ratio2)
})

test_that("inverse-Gamma log-likelihood ratios.", {
  skip_if_not_installed("invgamma")
  library(invgamma)

  x1 <- 2
  x2 <- 7
  shape <- 2
  rate <- 6

  ratio1 <- invGammaLogLikelihood(x1, shape, rate) - invGammaLogLikelihood(x2, shape, rate)
  ratio2 <- dinvgamma(x1, shape, rate, log = T) - dinvgamma(x2, shape, rate, log = T)

  expect_equal(ratio1, ratio2)
})

test_that("Wishart log-likelihood ratios.", {
  skip_if_not_installed("CholWishart")
  library(CholWishart)

  x1 <- matrix(c(3, 1, 1, 3), nrow = 2)
  x2 <- matrix(c(2, -0.7, -0.7, 2), nrow = 2)

  Psi <- matrix(c(1.23, 0.3, 0.3, 1.23), nrow = 2)

  nu <- 15
  P <- 2

  ratio1 <- wishartLogLikelihood(x1, Psi, nu, P) - wishartLogLikelihood(x2, Psi, nu, P)
  ratio2 <- dWishart(x1, nu, Psi, log = T) - dWishart(x2, nu, Psi, log = T)

  expect_equal(ratio1, ratio2)
})

test_that("Inverse-Wishart log-likelihood ratios.", {
  skip_if_not_installed("CholWishart")
  library(CholWishart)

  x1 <- matrix(c(3, 1, 1, 3), nrow = 2)
  x2 <- matrix(c(2, -0.7, -0.7, 2), nrow = 2)

  Psi <- matrix(c(1.23, 0.3, 0.3, 1.23), nrow = 2)

  nu <- 15
  P <- 2

  ratio1 <- invWishartLogLikelihood(x1, Psi, nu, P) - invWishartLogLikelihood(x2, Psi, nu, P)
  ratio2 <- dInvWishart(x1, nu, Psi, log = T) - dInvWishart(x2, nu, Psi, log = T)

  expect_equal(ratio1, ratio2)
})

test_that("LKJ prior sampler returns valid correlation matrices.", {
  # P is capped at 5 here: the rejection sampler's acceptance rate shrinks
  # combinatorially with P (see genericFunctions.cpp), so larger P is slow
  # and is not exercised in this fast unit test.
  for (P in 2:5) {
    for (eta in c(1, 2, 5)) {
      for (i in 1:20) {
        R <- sampleLKJCorrelationMatrix(P, eta)

        expect_equal(dim(R), c(P, P))
        expect_equal(R, t(R))
        expect_equal(diag(R), rep(1, P))
        expect_true(all(eigen(R, symmetric = TRUE, only.values = TRUE)$values > 0))
      }
    }
  }
})

test_that("LKJ prior sampler concentrates towards the identity as eta grows.", {
  set.seed(1)
  P <- 4
  n <- 300

  mean_abs_corr <- function(eta) {
    draws <- sapply(1:n, function(i) {
      R <- sampleLKJCorrelationMatrix(P, eta)
      mean(abs(R[upper.tri(R)]))
    })
    mean(draws)
  }

  # Higher eta shrinks correlations towards 0 (Lewandowski, Kurowicka & Joe,
  # 2009); check this holds in expectation across an increasing sequence.
  m1 <- mean_abs_corr(1)
  m5 <- mean_abs_corr(5)
  m20 <- mean_abs_corr(20)

  expect_true(m1 > m5)
  expect_true(m5 > m20)
})

test_that("LKJ prior sampler rejects eta < 1.", {
  expect_error(sampleLKJCorrelationMatrix(3, 0.5))
})
