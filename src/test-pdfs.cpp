# include <testthat.h>
# include "pdfs.h"
# include <RcppArmadillo.h>

bool compareDoubles(double x, double y, double epsilon = 1e-6) {
  double diff = x - y;
  return (diff < epsilon) && (-diff < epsilon);
}

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("Unit test for reduced PDF functions.") {
  
  test_that("gamma log likelihood") {
    expect_true(compareDoubles(gammaLogLikelihood(12, 21, 11), -74.2816827362277));
  }
  
  test_that("inverse gamma log likelihood") {
    expect_true(compareDoubles(invGammaLogLikelihood(21, 3, 17), -5.1811207088088));
  }
  
  arma::mat X, Psi;
  
  X = { {7, 2}, 
        {2, 7} };
  
  Psi = { {1.01, 0.50},
          {0.50, 1.01} };
  
  test_that("wishart log likelihood") {
    expect_true(compareDoubles(wishartLogLikelihood(X, Psi, 2, 2), -9.52418957709308));
  }
  
  test_that("inverse wishart log likelihood") {
    expect_true(compareDoubles(invWishartLogLikelihood(X, Psi, 2, 2), -9.39031021087776));
  }
  
}
