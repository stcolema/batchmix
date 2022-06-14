# 
# Unit tests for the log-likelihood functions implemented in C++.
# 
context("Likelihood comparison")
library(BatchMixtureModel)

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