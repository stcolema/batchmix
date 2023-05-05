
# include "genericFunctions.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

//' title Propose new non-negative value
//' description Propose new non-negative for sampling.
//' param x Current value to be proposed
//' param window The proposal window
//' return new double
double proposeNewNonNegativeValue(double x, double window, 
                                  bool use_log_norm,
                                  double tolerance
) {
  bool value_below_tolerance = false;
  double proposed_value = 0.0;
  if(use_log_norm) {
    proposed_value = std::exp(std::log(x) + randn() * window);
  } else {
    proposed_value = rGamma(x * window, window);
  }
  
  // If the value is too small (normally close to 0 or negative somehow)
  value_below_tolerance = (proposed_value < tolerance);
  if(value_below_tolerance) {
    proposed_value = proposeNewNonNegativeValue(x, window, use_log_norm, tolerance);
  }
  
  return proposed_value;
};

//' title The Inverse Gamma Distribution
//' description Random generation from the inverse Gamma distribution.
//' param shape Shape parameter.
//' param rate Rate parameter.
//' return Sample from invGamma(shape, rate).
double rInvGamma(double shape, double rate) {
  double x = arma::randg( distr_param(shape, 1.0 / rate) );
  return (1 / x);
};

//' title The Inverse Gamma Distribution
//' description Random generation from the inverse Gamma distribution.
//' param N Number of samples to draw.
//' param shape Shape parameter.
//' param rate Rate parameter.
//' return Sample from invGamma(shape, rate).
arma::vec rInvGamma(uword N, double shape, double rate) {
  vec x = arma::randg(N, distr_param(shape, 1.0 / rate) );
  return (1 / x);
};

//' title The Gamma Distribution
//' description Random generation from the Gamma distribution.
//' param shape Shape parameter.
//' param rate Rate parameter.
//' return Sample from Gamma(shape, rate).
double rGamma(double shape, double rate) {
  return arma::randg( distr_param(shape, 1.0 / rate) );
};

//' title The Gamma Distribution
//' description Random generation from the Gamma distribution.
//' param N Number of samples to draw.
//' param shape Shape parameter.
//' param rate Rate parameter.
//' return N samples from Gamma(shape, rate).
arma::vec rGamma(uword N, double shape, double rate) {
  return arma::randg(N, distr_param(shape, 1.0 / rate) );
};

//' title The Beta Distribution
//' description Random generation from the Beta distribution.
//' See https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions.
//' Samples from a Beta distribution based using two independent gamma
//' distributions.
//' param a Shape parameter.
//' param b Shape parameter.
//' return Sample from Beta(a, b).
double rBeta(double a, double b) { // double theta = 1.0) {
  double X = arma::randg( arma::distr_param(a, 1.0) );
  double Y = arma::randg( arma::distr_param(b, 1.0) );
  double beta = X / (double)(X + Y);
  return(beta);
};

//' title The Beta Distribution
//' description Random generation from the Beta distribution.
//' See https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions.
//' Samples from a Beta distribution based using two independent gamma
//' distributions.
//' param n The number of samples to draw.
//' param a Shape parameter.
//' param b Shape parameter.
//' return Sample from Beta(a, b).
arma::vec rBeta(arma::uword n, double a, double b) {
  arma::vec X = arma::randg(n, arma::distr_param(a, 1.0) );
  arma::vec Y = arma::randg(n, arma::distr_param(b, 1.0) );
  arma::vec beta = X / (X + Y);
  return(beta);
};

//' title Metropolis acceptance step
//' description Given a probaility, randomly accepts by sampling from a uniform 
//' distribution.
//' param acceptance_prob Double between 0 and 1.
//' return Boolean indicating acceptance.
bool metropolisAcceptanceStep(double acceptance_prob) {
  double u = arma::randu();
  return (u < acceptance_prob);
};

//' title Sample mean
//' description calculate the sample mean of a matrix X.
//' param X Matrix
//' return Vector of the column means of X.
vec sampleMean(arma::mat X) {
  mat mu_t = mean(X);
  return mu_t.row(0).t();
};

//' title Calculate sample covariance
//' description Returns the unnormalised sample covariance. Required as
//' arma::cov() does not work for singletons.
//' param data Data in matrix format
//' param sample_mean Sample mean for data
//' param n The number of samples in data
//' param n_col The number of columns in data
//' return One of the parameters required to calculate the posterior of the
//'  Multivariate normal with uknown mean and covariance (the unnormalised
//'  sample covariance).
arma::mat calcSampleCov(arma::mat data,
                        arma::vec sample_mean,
                        arma::uword N,
                        arma::uword P
) {
  
  mat sample_covariance = zeros<mat>(P, P);
  
  // If n > 0 (as this would crash for empty clusters), and for n = 1 the
  // sample covariance is 0
  if(N > 1){
    data.each_row() -= sample_mean.t();
    sample_covariance = data.t() * data;
  }
  return sample_covariance;
};
