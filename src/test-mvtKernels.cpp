// mvnSampler.cpp
// =============================================================================
//
// included dependencies
# include <RcppArmadillo.h>
# include <testthat.h>
# include "mvtSampler.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// mvnSampler unit tests

bool compareDoubles3(double x, double y, double epsilon = 1e-6) {
  double diff = x - y;
  return (diff < epsilon) && (-diff < epsilon);
}

context("Unit test for MVT posterior kernels.") {
  
  uword K = 2, B = 3;
  double
    mu_proposal_window = 0.5,
      cov_proposal_window = 200,
      m_proposal_window = 0.4,
      S_proposal_window = 100,
      t_df_proposal_window = 35,
      m_scale = 0.01,
      rho = 3.0,
      theta = 1.0;
  
  uvec labels(10), batch_vec(10);
  vec concentration(K);

  labels = {0, 1, 1, 0, 1, 1, 1, 0, 0, 1};
  batch_vec = {0, 1, 0, 1, 0, 1, 1, 2, 2, 2};
  concentration = {1.0, 1.0};

  arma::mat Y = { 7.2,
    3.1,
    2.2,
    9.8,
    2.3,
    3.8,
    3.3,
    5.2,
    6.8,
    1.3
  }, X = Y.t();
  
  mvtSampler toy_sampler(
    K,
    B,
    mu_proposal_window,
    cov_proposal_window,
    m_proposal_window,
    S_proposal_window,
    t_df_proposal_window,
    labels,
    batch_vec,
    concentration,
    X,
    m_scale,
    rho,
    theta
  );
  
  double val1 = 0.0,
    val2 = 0.0,
    val3 = 0.0,
    val4 = 0.0,
    val5 = 0.0,
    val6 = 0.0,
    val7 = 0.0,
    val8 = 0.0,
    val9 = 0.0,
    val10 = 0.0,
    val11 = 0.0,
    val12 = 0.0;
  
  vec mu_0 = {3.0}, mu_1 = {7.2};
  mat m = { {0.0, 1.0, -1.0} },
    mean_sum_0 = m + 3.0,
    mean_sum_1 = m + 7.2;
  
  // Initialise some mean vectors
  toy_sampler.mu.col(0) = {3.0};
  toy_sampler.mu.col(1) = {7.2};
  
  toy_sampler.m.col(0) = {0.0};
  toy_sampler.m.col(1) = {1.0};
  toy_sampler.m.col(2) = {-1.0};
  
  // Initialise some covariance matrices;
  toy_sampler.cov.slice(0) = {1.0};
  toy_sampler.cov.slice(1) = {0.7};
  
  toy_sampler.S.col(0) = {1.2};
  toy_sampler.S.col(1) = {1.3};
  toy_sampler.S.col(2) = {1.5};
  
  toy_sampler.t_df(0) = 13;
  toy_sampler.t_df(1) = 49;
  
  // Matrix combinations
  toy_sampler.matrixCombinations();
  
  val1 = toy_sampler.muLogKernel(0, mu_0, mean_sum_0);
  val2 = toy_sampler.muLogKernel(1, mu_1, mean_sum_1);

  test_that("mu log posterior kernel") {
    expect_true(compareDoubles3(val1, -22.017362, 1e-5));
    expect_true(compareDoubles3(val2, -63.961374, 1e-5));
  }
  
  val3 = toy_sampler.mLogKernel(0, m.col(0), toy_sampler.mean_sum.cols(0, 1) );
  val4 = toy_sampler.mLogKernel(1, m.col(1), toy_sampler.mean_sum.cols(2, 3) );
  val5 = toy_sampler.mLogKernel(2, m.col(2), toy_sampler.mean_sum.cols(4, 5) );
  
  test_that("m log posterior kernel") {
    expect_true(compareDoubles3(val3, -8.802964, 1e-5));
    expect_true(compareDoubles3(val4, -38.648574, 1e-5));
    expect_true(compareDoubles3(val5, -20.210158, 1e-5));
  }
  
  val6 = toy_sampler.covLogKernel(0, 
    toy_sampler.cov.slice(0), 
    toy_sampler.cov_log_det(0),
    toy_sampler.cov_inv.slice(0),
    toy_sampler.cov_comb_log_det.row(0).t(),
    toy_sampler.cov_comb_inv.slices(0, 2)
  );
  
  val7 = toy_sampler.covLogKernel(1, 
    toy_sampler.cov.slice(1), 
    toy_sampler.cov_log_det(1),
    toy_sampler.cov_inv.slice(1),
    toy_sampler.cov_comb_log_det.row(1).t(),
    toy_sampler.cov_comb_inv.slices(3, 5)
  );
  
  test_that("cov log posterior kernel") {
    expect_true(compareDoubles3(val6, -22.925973, 1e-5));
    expect_true(compareDoubles3(val7, -64.189365, 1e-5));
  }
  
  val8 = toy_sampler.sLogKernel(0,
    toy_sampler.S.col(0),
    toy_sampler.cov_comb_log_det.col(0),
    toy_sampler.cov_comb_inv.slices(0, 3)
  );
  
  val9 = toy_sampler.sLogKernel(1,
    toy_sampler.S.col(1),
    toy_sampler.cov_comb_log_det.col(1),
    toy_sampler.cov_comb_inv.slices(1, 4)
  );
  
  val10 = toy_sampler.sLogKernel(2,
    toy_sampler.S.col(2),
    toy_sampler.cov_comb_log_det.col(2),
    toy_sampler.cov_comb_inv.slices(2, 5)
  );
  
  test_that("s log posterior kernel") {
    expect_true(compareDoubles3(val8, -36.077722, 1e-5));
    expect_true(compareDoubles3(val9, -34.990054, 1e-5));
    expect_true(compareDoubles3(val10, -25.083692, 1e-5));
  }
  
  val11 = toy_sampler.dfLogKernel(0, toy_sampler.t_df(0), toy_sampler.pdf_coef(0));
  val12 = toy_sampler.dfLogKernel(1, toy_sampler.t_df(1), toy_sampler.pdf_coef(1));
  
  test_that("df log posterior kernel") {
    expect_true(compareDoubles3(val11, -20.553550, 1e-5));
    expect_true(compareDoubles3(val12, -62.930631, 1e-5));
  }
}
