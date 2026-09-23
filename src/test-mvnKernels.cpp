// mvnSampler.cpp
// =============================================================================
//
// included dependencies
# include <RcppArmadillo.h>
# include <testthat.h>
# include "mvnSampler.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// mvnSampler unit tests

bool compareDoubles2(double x, double y, double epsilon = 1e-6) {
  double diff = x - y;
  return (diff < epsilon) && (-diff < epsilon);
}

context("Unit test for MVN posterior kernels.") {

  bool sample_m_scale = false;
  uword K = 2, B = 3;
  double
    mu_proposal_window = 0.5,
      cov_proposal_window = 200,
      m_proposal_window = 0.4,
      S_proposal_window = 100,
      m_scale = 0.01,
      rho = 3.0,
      theta = 1.0;

  uvec labels(10), batch_vec(10), fixed(10, fill::zeros);
  vec concentration(K);
  // mat X(10, 1);

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

  mvnSampler toy_sampler(
      K,
      B,
      mu_proposal_window,
      cov_proposal_window,
      m_proposal_window,
      S_proposal_window,
      labels,
      batch_vec,
      concentration,
      X,
      fixed,
      m_scale,
      rho,
      theta,
      sample_m_scale
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
    val10 = 0.0;

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

  // For the covariance combination
  toy_sampler.matrixCombinations();

  val1 = toy_sampler.muLogKernel(0, mu_0, mean_sum_0);
  val2 = toy_sampler.muLogKernel(1, mu_1, mean_sum_1);
  
  test_that("mu log posterior kernel") {
    expect_true(compareDoubles2(val1, -32.02085, 1e-5));
    expect_true(compareDoubles2(val2, -78.48748, 1e-5));
  }

  val3 = toy_sampler.mLogKernel(0, m.col(0), toy_sampler.mean_sum.cols(0, 1) );
  val4 = toy_sampler.mLogKernel(1, m.col(1), toy_sampler.mean_sum.cols(2, 3) );
  val5 = toy_sampler.mLogKernel(2, m.col(2), toy_sampler.mean_sum.cols(4, 5) );

  test_that("m log posterior kernel") {
    expect_true(compareDoubles2(val3, -10.91562, 1e-5));
    expect_true(compareDoubles2(val4, -54.2134, 1e-5));
    expect_true(compareDoubles2(val5, -22.39516, 1e-5));
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
    expect_true(compareDoubles2(val6, -32.92946, 1e-5));
    expect_true(compareDoubles2(val7, -78.71547, 1e-5));
  }
  
  val8 = toy_sampler.sLogKernel(0,
    toy_sampler.S.col(0),
    toy_sampler.cov_comb_log_det.col(0),
    toy_sampler.cov_comb_inv.slices(toy_sampler.KB_inds + 0)
  );

  val9 = toy_sampler.sLogKernel(1,
    toy_sampler.S.col(1),
    toy_sampler.cov_comb_log_det.col(1),
    toy_sampler.cov_comb_inv.slices(toy_sampler.KB_inds + 1)
  );

  val10 = toy_sampler.sLogKernel(2,
    toy_sampler.S.col(2),
    toy_sampler.cov_comb_log_det.col(2),
    toy_sampler.cov_comb_inv.slices(toy_sampler.KB_inds + 2)
  );
  
  test_that("s log posterior kernel") {
    expect_true(compareDoubles2(val8, -35.00168, 1e-5));
    expect_true(compareDoubles2(val9, -49.56650, 1e-5));
    expect_true(compareDoubles2(val10, -22.18394, 1e-5));
  }


}

// =============================================================================
// Regression: sampleMScalePosterior()/sampleMPrior() used the wrong scale
// conventions. sampleMScalePosterior() computed the InvGamma posterior rate
// as sum(m^2)/(4*delta_2) instead of the correct sum(m^2)/(2*delta_2) (an
// erroneous extra factor of 0.5); sampleMPrior() drew m ~ N(mean, precision)
// instead of N(mean, 1/precision) (batch_shift_prior_precision is a
// precision, so the standard deviation is its inverse square root, not the
// precision itself). Both are law-of-large-numbers checks: fix
// delta_2/lambda_2/m to known values, redraw many times, and compare the
// empirical mean/variance against the closed-form InvGamma/Normal moments
// the correct formula implies - tight enough that the old bugs (rate off by
// a factor of 2; variance off by a factor of precision^3) would fail by a
// wide margin, loose enough to tolerate ordinary Monte Carlo noise.
context("Regression: sampleMScalePosterior/sampleMPrior scale conventions.") {

  bool sample_m_scale = true;
  uword K = 2, B = 3;
  double
    mu_proposal_window = 0.5,
      cov_proposal_window = 200,
      m_proposal_window = 0.4,
      S_proposal_window = 100,
      m_scale = 0.01,
      rho = 3.0,
      theta = 1.0;

  uvec labels(10), batch_vec(10), fixed(10, fill::zeros);
  vec concentration(K);

  labels = {0, 1, 1, 0, 1, 1, 1, 0, 0, 1};
  batch_vec = {0, 1, 0, 1, 0, 1, 1, 2, 2, 2};
  concentration = {1.0, 1.0};

  arma::mat Y = { 7.2, 3.1, 2.2, 9.8, 2.3, 3.8, 3.3, 5.2, 6.8, 1.3 }, X = Y.t();

  mvnSampler scale_sampler(
      K, B, mu_proposal_window, cov_proposal_window, m_proposal_window,
      S_proposal_window, labels, batch_vec, concentration, X, fixed,
      m_scale, rho, theta, sample_m_scale
  );

  // sampleMScalePosterior(): a_pos = a + 0.5*P*B was already correct; the
  // bug was in b_pos, which should be b + sum(m^2)/(2*delta_2).
  scale_sampler.delta_2 = 1.0;
  scale_sampler.m = { {0.0, 1.0, -1.0} }; // P = 1, B = 3; sum(m^2) = 2.0

  double a_pos_expected = scale_sampler.a + 0.5 * (double) scale_sampler.P * (double) B;
  double b_pos_expected = scale_sampler.b + 2.0 / (2.0 * scale_sampler.delta_2);
  double expected_lambda2_mean = b_pos_expected / (a_pos_expected - 1.0);

  uword n_draws = 40000;
  double lambda2_sum = 0.0;
  for (uword i = 0; i < n_draws; i++) {
    scale_sampler.sampleMScalePosterior();
    lambda2_sum += scale_sampler.lambda_2;
  }
  double lambda2_mean = lambda2_sum / (double) n_draws;

  test_that("sampleMScalePosterior draws lambda_2 from the correctly-scaled InvGamma posterior") {
    // Loose (5%) relative tolerance for Monte Carlo noise; the old
    // (halved-rate) bug would put the true posterior mean roughly a factor
    // of 2 away from this, far outside a 5% band.
    expect_true(compareDoubles2(lambda2_mean, expected_lambda2_mean, 0.05 * expected_lambda2_mean));
  }

  // sampleMPrior(): m ~ N(mean, 1/precision), so its standard deviation is
  // sqrt(1/precision), not precision itself.
  scale_sampler.batch_shift_prior_precision = 0.25; // Var should be 1/0.25 = 4.0
  scale_sampler.batch_shift_prior_mean = 0.0;

  double m_sum = 0.0, m_sq_sum = 0.0;
  uword n_m_draws = 40000;
  for (uword i = 0; i < n_m_draws; i++) {
    scale_sampler.sampleMPrior();
    for (uword bb = 0; bb < B; bb++) {
      double v = scale_sampler.m(0, bb);
      m_sum += v;
      m_sq_sum += v * v;
    }
  }
  double n_total = (double) (n_m_draws * B);
  double m_mean = m_sum / n_total;
  double m_var = m_sq_sum / n_total - m_mean * m_mean;

  test_that("sampleMPrior draws m with variance 1/precision, not precision^2") {
    // True variance is 4.0; the old bug would give precision^2 = 0.0625 -
    // two orders of magnitude away, so a generous absolute tolerance is
    // both tight enough to catch the bug and loose enough for Monte Carlo
    // noise on a variance estimate from 3 * 40000 draws.
    expect_true(compareDoubles2(m_var, 4.0, 0.6));
    expect_true(compareDoubles2(m_mean, 0.0, 0.1));
  }
}
