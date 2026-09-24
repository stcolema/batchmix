
# include "genericFunctions.h"
#include <cmath>
#include <ctime>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::mat imputeColumnMeans(arma::mat X) {
  for(uword p = 0; p < X.n_cols; p++) {
    vec col_p = X.col(p);
    uvec observed = find_finite(col_p);
    double col_mean = observed.n_elem > 0 ? mean(col_p.elem(observed)) : 0.0;
    for(uword n = 0; n < X.n_rows; n++) {
      if(!std::isfinite(X(n, p))) {
        X(n, p) = col_mean;
      }
    }
  }
  return X;
};

//' @title Propose new non-negative value
//' @description Propose new non-negative for sampling.
//' @param x Current value to be proposed
//' @param window The proposal window
//' @return new double
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

//' @title The Inverse Gamma Distribution
//' @description Random generation from the inverse Gamma distribution.
//' @param shape Shape parameter.
//' @param rate Rate parameter.
//' @return Sample from invGamma(shape, rate).
double rInvGamma(double shape, double rate) {
  double x = arma::randg( distr_param(shape, 1.0 / rate) );
  return (1 / x);
};

//' @title The Inverse Gamma Distribution
//' @description Random generation from the inverse Gamma distribution.
//' @param N Number of samples to draw.
//' @param shape Shape parameter.
//' @param rate Rate parameter.
//' @return Sample from invGamma(shape, rate).
arma::vec rInvGamma(uword N, double shape, double rate) {
  vec x = arma::randg(N, distr_param(shape, 1.0 / rate) );
  return (1 / x);
};


//' @title The Gamma Distribution
//' @description Random generation from the Gamma distribution.
//' @param shape Shape parameter.
//' @param rate Rate parameter.
//' @return Sample from Gamma(shape, rate).
double rGamma(double shape, double rate) {
  return arma::randg( distr_param(shape, 1.0 / rate) );
};

//' @title The Gamma Distribution
//' @description Random generation from the Gamma distribution.
//' @param N Number of samples to draw.
//' @param shape Shape parameter.
//' @param rate Rate parameter.
//' @return N samples from Gamma(shape, rate).
arma::vec rGamma(uword N, double shape, double rate) {
  return arma::randg(N, distr_param(shape, 1.0 / rate) );
};

//' @title The Beta Distribution
//' @description Random generation from the Beta distribution.
//' See https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions.
//' Samples from a Beta distribution based using two independent gamma
//' distributions.
//' @param a Shape parameter.
//' @param b Shape parameter.
//' @return Sample from Beta(a, b).
double rBeta(double a, double b) { // double theta = 1.0) {
  double X = arma::randg( arma::distr_param(a, 1.0) );
  double Y = arma::randg( arma::distr_param(b, 1.0) );
  double beta = X / (double)(X + Y);
  return(beta);
};

//' @title The Beta Distribution
//' @description Random generation from the Beta distribution.
//' See https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions.
//' Samples from a Beta distribution based using two independent gamma
//' distributions.
//' @param n The number of samples to draw.
//' @param a Shape parameter.
//' @param b Shape parameter.
//' @return Sample from Beta(a, b).
arma::vec rBeta(arma::uword n, double a, double b) {
  arma::vec X = arma::randg(n, arma::distr_param(a, 1.0) );
  arma::vec Y = arma::randg(n, arma::distr_param(b, 1.0) );
  arma::vec beta = X / (X + Y);
  return(beta);
};

//' @title The Log-Normal Distribution
 //' description Random generation from the log-Normal distribution.
 //' param mu mean parameter.
 //' param sd standard deviation parameter.
 //' return Sample from log-Normal(mu, sd^2).
 double rLogNormal(double mu, double sd) {
   return log(arma::randn<double>( distr_param(mu, sd) ));
 };
 
 //' title The Log-Normal Distribution
 //' description Random generation from the log-Normal distribution.
 //' param N positive integer - the number of samples drawn.
 //' param mu mean parameter.
 //' param sd standard deviation parameter.
 //' return N samples from log-Normal(mu, sd^2)
 arma::vec rLogNormal(arma::uword N, double mu, double sd) {
   return arma::log(arma::randn<arma::vec>( N, distr_param(mu,sd) ));
 };

//' @title Metropolis acceptance step
//' @description Given a probaility, randomly accepts by sampling from a uniform 
//' distribution.
//' @param acceptance_prob Double between 0 and 1.
//' @return Boolean indicating acceptance.
bool metropolisAcceptanceStep(double acceptance_prob) {
  double u = arma::randu();
  return (u < acceptance_prob);
};

//' @title Accept proposal
//' @description Determines if a proposal is accepted given a log ratio of scores 
//' for the proposed and original values.
//' @param proposed_model_score Score in the posterior kernel for the proposed 
//' parameter value
//' @param current_model_score Score in the posterior kernel for the current 
//' parameter value
//' @return Boolean indicating acceptance.
bool acceptProposal(double proposed_model_score, double current_model_score) {
  double u = randu(), acceptance_prob = 0.0;
  acceptance_prob = std::min(
    1.0, 
    std::exp(proposed_model_score - current_model_score)
  );
  return u < acceptance_prob;
};

//' @title Sample mean
//' @description calculate the sample mean of a matrix X.
//' @param X Matrix
//' @return Vector of the column means of X.
vec sampleMean(arma::mat X) {
  mat mu_t = mean(X);
  return mu_t.row(0).t();
};

//' @title Calculate sample covariance
//' @description Returns the unnormalised sample covariance. Required as
//' arma::cov() does not work for singletons.
//' @param data Data in matrix format
//' @param sample_mean Sample mean for data
//' @param n The number of samples in data
//' @param n_col The number of columns in data
//' @return One of the parameters required to calculate the posterior of the
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

//' @title The LKJ distribution
//' @description Random generation of a correlation matrix from LKJ(eta) by
//' rejection sampling; see the header for the derivation.
//' @param P Dimension of the correlation matrix.
//' @param eta Concentration parameter, eta >= 1.
//' @return A P x P correlation matrix sampled from LKJ(eta).
// [[Rcpp::export]]
arma::mat sampleLKJCorrelationMatrix(arma::uword P, double eta) {

  if(eta < 1.0) {
    Rcpp::stop("sampleLKJCorrelationMatrix: eta < 1 is not supported (the LKJ density is unbounded near singular correlation matrices for eta < 1).");
  }

  mat R = eye<mat>(P, P);
  if(P < 2) {
    return R;
  }

  mat candidate(P, P), L(P, P);
  bool is_pd = false, accepted = false;
  double det_r = 0.0, accept_prob = 0.0;
  uword n_attempts = 0;
  const uword max_attempts = 2000000;

  while(!accepted) {
    n_attempts++;
    if(n_attempts > max_attempts) {
      Rcpp::stop("sampleLKJCorrelationMatrix: exceeded " + std::to_string(max_attempts) + " rejection-sampling attempts for P = " + std::to_string(P) + ", eta = " + std::to_string(eta) + ". Rejection sampling from a uniform box becomes impractical for larger P (the fraction of the box that is positive definite shrinks combinatorially); this method is only suitable for small P.");
    }
    candidate = eye<mat>(P, P);
    for(uword i = 0; i < P; i++) {
      for(uword j = i + 1; j < P; j++) {
        double r_ij = 2.0 * randu() - 1.0;
        candidate(i, j) = r_ij;
        candidate(j, i) = r_ij;
      }
    }

    // Reject candidates that are not valid (positive definite) correlation
    // matrices; the non-throwing form of chol() returns false rather than
    // raising an exception on failure.
    is_pd = arma::chol(L, candidate);
    if(!is_pd) {
      continue;
    }

    det_r = arma::det(candidate);
    if(det_r <= 0.0) {
      continue;
    }

    if(eta == 1.0) {
      accepted = true;
    } else {
      accept_prob = std::pow(det_r, eta - 1.0);
      accepted = (randu() < accept_prob);
    }

    if(accepted) {
      R = candidate;
    }
  }

  return R;
};

//' @title The truncated Normal distribution (right-tail helper)
//' @description Robert (1995, "Simulation of truncated normal variables",
//' Statistics and Computing 5(2)) exponential-tilting rejection sampler
//' for a standard Normal truncated to (alpha, infinity), used when alpha
//' is far enough into the tail that inverse-CDF sampling loses precision
//' (the CDF saturates to 1 in floating point). Returns a draw of the
//' STANDARDISED variable, i.e. already on the (lower - mean)/sd scale.
double rTruncNormRightTailStd(double alpha) {
  double a_star = 0.5 * (alpha + std::sqrt(alpha * alpha + 4.0));
  double z = 0.0, rho = 0.0;
  bool accepted = false;
  while(!accepted) {
    z = alpha - std::log(randu()) / a_star;
    rho = std::exp(-0.5 * std::pow(z - a_star, 2.0));
    accepted = (randu() <= rho);
  }
  return z;
};

//' @title The truncated Normal distribution
//' @description Random generation from a truncated Normal. Uses inverse
//' CDF sampling in the regime where that is numerically reliable, and
//' falls back to Robert's (1995) exponential-tilting rejection sampler
//' (see rTruncNormRightTailStd()) for one-sided truncation far into a
//' tail, where inverse CDF sampling would otherwise saturate to the
//' truncation boundary itself rather than a proper draw. Two-sided
//' truncation with both bounds simultaneously far into the same tail is
//' not specially handled and falls back to inverse CDF, per the header.
//' @param mean Mean of the untruncated Normal distribution.
//' @param sd Standard deviation of the untruncated Normal distribution.
//' @param lower Lower truncation bound (-arma::datum::inf for none).
//' @param upper Upper truncation bound (arma::datum::inf for none).
//' @return A draw from Normal(mean, sd^2) truncated to (lower, upper).
double rTruncNorm(double mean, double sd, double lower, double upper) {

  const double tail_threshold = 5.0;
  bool lower_is_inf = std::isinf(lower);
  bool upper_is_inf = std::isinf(upper);

  if(upper_is_inf && !lower_is_inf) {
    double alpha = (lower - mean) / sd;
    if(alpha > tail_threshold) {
      return mean + sd * rTruncNormRightTailStd(alpha);
    }
  }

  if(lower_is_inf && !upper_is_inf) {
    double beta = (upper - mean) / sd;
    if(beta < -tail_threshold) {
      // Reflect about the mean: X ~ TruncNorm(mean, sd, -inf, upper) has
      // the same distribution as 2*mean - Y for Y ~ TruncNorm(mean, sd,
      // 2*mean - upper, inf), which is the right-tail case above.
      double alpha = -beta;
      return mean - sd * rTruncNormRightTailStd(alpha);
    }
  }

  double p_lower = lower_is_inf ? 0.0 : R::pnorm(lower, mean, sd, 1, 0);
  double p_upper = upper_is_inf ? 1.0 : R::pnorm(upper, mean, sd, 1, 0);

  // Guard against a degenerate (zero-width, in floating point) interval;
  // this can still occur for two-sided far-tail truncation, which is not
  // covered by the rejection sampler above.
  if(p_upper <= p_lower) {
    return lower_is_inf ? upper : lower;
  }

  double u = p_lower + randu() * (p_upper - p_lower);

  // Keep u strictly inside (0, 1) so qnorm does not return +/-Inf.
  u = std::min(std::max(u, 1e-12), 1.0 - 1e-12);

  return R::qnorm(u, mean, sd, 1, 0);
};

//' @title Squared-exponential covariance kernel
//' @description Builds a Gaussian process covariance matrix; see header.
//' @param x Vector of 1-D locations.
//' @param tau2 Marginal variance.
//' @param length_scale Correlation length scale.
//' @param jitter Diagonal jitter for numerical stability.
//' @return The covariance matrix.
// [[Rcpp::export]]
arma::mat squaredExponentialKernel(arma::vec x, double tau2, double length_scale, double jitter) {

  uword n = x.n_elem;
  mat K(n, n);
  double d = 0.0;

  for(uword i = 0; i < n; i++) {
    for(uword j = 0; j < n; j++) {
      d = x(i) - x(j);
      K(i, j) = tau2 * std::exp(-(d * d) / (2.0 * length_scale * length_scale));
    }
  }
  K.diag() += jitter;

  return K;
};

//' @title Matern-3/2 covariance kernel
//' @description Builds a Gaussian process covariance matrix; see header
//' for the full derivation, references and the conditioning rationale for
//' preferring this over \code{squaredExponentialKernel()}.
//' @param x Vector of 1-D locations.
//' @param tau2 Marginal variance.
//' @param length_scale Correlation length scale.
//' @param jitter Diagonal jitter for numerical stability.
//' @return The covariance matrix.
// [[Rcpp::export]]
arma::mat maternKernel32(arma::vec x, double tau2, double length_scale, double jitter) {

  uword n = x.n_elem;
  mat K(n, n);
  double d = 0.0, scaled_d = 0.0;
  const double sqrt3 = std::sqrt(3.0);

  for(uword i = 0; i < n; i++) {
    for(uword j = 0; j < n; j++) {
      d = std::abs(x(i) - x(j));
      scaled_d = sqrt3 * d / length_scale;
      K(i, j) = tau2 * (1.0 + scaled_d) * std::exp(-scaled_d);
    }
  }
  K.diag() += jitter;

  return K;
};

//' @title Multinomial-logit Gaussian process log-kernel
//' @description The unnormalised log-posterior-kernel for one ALR
//' coordinate of batch-dependent multinomial weights under a GP prior
//' over the batch index; see header for the full derivation and
//' references.
//' @param eta B-vector, this ALR coordinate for each batch.
//' @param eta_other_sum B-vector, the softmax normalising contribution of
//' every other non-pivot category, held fixed this step.
//' @param class_counts_j B-vector, per-batch counts in this category.
//' @param class_counts_total B-vector, per-batch total item counts.
//' @return The log-likelihood value for eta.
// [[Rcpp::export]]
double multinomialLogitLogLik(
  arma::vec eta,
  arma::vec eta_other_sum,
  arma::vec class_counts_j,
  arma::vec class_counts_total
) {
  double log_lik = 0.0, D_b = 0.0;

  for(uword b = 0; b < eta.n_elem; b++) {
    D_b = 1.0 + std::exp(eta(b)) + eta_other_sum(b);
    log_lik += class_counts_j(b) * eta(b) - class_counts_total(b) * std::log(D_b);
  }

  return log_lik;
};

//' @param eta B-vector, this ALR coordinate for each batch.
//' @param eta_other_sum B-vector, the softmax normalising contribution of
//' every other non-pivot category, held fixed this step.
//' @param class_counts_j B-vector, per-batch counts in this category.
//' @param class_counts_total B-vector, per-batch total item counts.
//' @param gp_chol The LOWER-triangular Cholesky factor of the GP
//' covariance matrix for this coordinate; see header for why a triangular
//' solve against this, rather than forming/using the covariance's inverse
//' directly, is the numerically-preferred modern approach.
//' @param beta This coordinate's estimated GP intercept; see header.
//' @return The unnormalised log-posterior-kernel value for eta.
// [[Rcpp::export]]
double multinomialLogitGPLogKernel(
  arma::vec eta,
  arma::vec eta_other_sum,
  arma::vec class_counts_j,
  arma::vec class_counts_total,
  arma::mat gp_chol,
  double beta
) {

  double log_lik = multinomialLogitLogLik(eta, eta_other_sum, class_counts_j, class_counts_total);

  // Centred on beta, not 0 - see header for why a fixed zero-mean GP is a
  // needlessly restrictive special case of the model this is meant to
  // implement. v' Sigma^-1 v = ||L^-1 v||^2 exactly (Sigma = L L'), via a
  // forward triangular solve - never Sigma^-1 itself; see header.
  vec centred = eta - beta;
  vec whitened = arma::solve(arma::trimatl(gp_chol), centred);
  double log_prior = -0.5 * arma::dot(whitened, whitened);

  return log_lik + log_prior;
};

// Not Rcpp::export'd: reference (call-by-mutation) output parameters are
// not something Rcpp Attributes turns into a sensible R-callable wrapper -
// this is a C++-internal helper, called only from sampler.cpp.
void invGammaQuantileMatch(
  double lower, double upper, double tail_prob,
  double& shape_out, double& rate_out
) {
  double target_ratio = upper / lower;

  // ratio_for_shape(a) = qgamma(1 - tail_prob, a) / qgamma(tail_prob, a)
  // (both at Gamma(a, rate = 1), so the ratio is independent of the
  // Inverse-Gamma rate we are solving for) is monotonically DEcreasing in
  // a - verified numerically (a = 0.1 gives ~2.6e20; a = 100 gives ~1.6),
  // so a plain bisection on a is sufficient and robust, with no derivative
  // needed and no risk of divergence.
  double a_lo = 1e-2, a_hi = 1e3;
  for (int iter = 0; iter < 100; iter++) {
    double a_mid = 0.5 * (a_lo + a_hi);
    double q_hi = R::qgamma(1.0 - tail_prob, a_mid, 1.0, 1, 0);
    double q_lo = R::qgamma(tail_prob, a_mid, 1.0, 1, 0);
    double ratio_mid = q_hi / q_lo;
    if (ratio_mid > target_ratio) {
      a_lo = a_mid; // ratio still too big -> need a larger shape to shrink it
    } else {
      a_hi = a_mid;
    }
  }
  double shape = 0.5 * (a_lo + a_hi);

  // Given shape, back out rate from matching the lower quantile exactly:
  // lower = 1 / qgamma(1 - tail_prob, shape, rate = rate_out), and
  // qgamma(p, shape, rate = r) = qgamma(p, shape, scale = 1) / r, so
  // rate_out = lower * qgamma(1 - tail_prob, shape, scale = 1).
  double q_hi_final = R::qgamma(1.0 - tail_prob, shape, 1.0, 1, 0);
  double rate = lower * q_hi_final;

  shape_out = shape;
  rate_out = rate;
};

//' @title Build a correlation-matrix Cholesky factor from partial
//' correlations
//' @description See header for the construction and its provenance.
//' @param Z A P x P matrix; only strictly-lower-triangular entries used.
//' @param P The dimension.
//' @return The P x P lower-triangular Cholesky factor L.
// [[Rcpp::export]]
arma::mat buildCorrelationCholeskyFromZ(arma::mat Z, arma::uword P) {

  mat L = zeros<mat>(P, P);
  L(0, 0) = 1.0;

  for(uword i = 1; i < P; i++) {
    double running_sum = 0.0;
    for(uword j = 0; j < i; j++) {
      double remaining = std::sqrt(std::max(0.0, 1.0 - running_sum));
      L(i, j) = Z(i, j) * remaining;
      running_sum += L(i, j) * L(i, j);
    }
    L(i, i) = std::sqrt(std::max(0.0, 1.0 - running_sum));
  }

  return L;
};

//' @title Invert buildCorrelationCholeskyFromZ()
//' @description Recovers the partial correlations Z from a valid
//' correlation-matrix Cholesky factor L.
//' @param L The P x P Cholesky factor.
//' @param P The dimension.
//' @return The P x P matrix Z.
// [[Rcpp::export]]
arma::mat choleskyToPartialCorrelations(arma::mat L, arma::uword P) {

  mat Z = zeros<mat>(P, P);

  for(uword i = 1; i < P; i++) {
    double running_sum = 0.0;
    for(uword j = 0; j < i; j++) {
      double remaining = std::sqrt(std::max(1e-300, 1.0 - running_sum));
      Z(i, j) = L(i, j) / remaining;
      running_sum += L(i, j) * L(i, j);
    }
  }

  return Z;
};

//' @title Jacobian of the partial-correlation to correlation-matrix map
//' @description log|dR/dZ|; see header for verification against
//' finite-difference Jacobians.
//' @param Z A P x P matrix; only strictly-lower-triangular entries used.
//' @param P The dimension.
//' @return The log-Jacobian determinant.
// [[Rcpp::export]]
double logJacobianZToR(arma::mat Z, arma::uword P) {

  double log_jac = 0.0;
  mat L = zeros<mat>(P, P);
  L(0, 0) = 1.0;

  // The |dL/dZ| part: the per-row product of "remaining L2-norm budget"
  // factors used to build each row of L (see buildCorrelationCholeskyFromZ).
  for(uword i = 1; i < P; i++) {
    double running_sum = 0.0;
    for(uword j = 0; j < i; j++) {
      double remaining = std::sqrt(std::max(1e-300, 1.0 - running_sum));
      log_jac += std::log(remaining);
      L(i, j) = Z(i, j) * remaining;
      running_sum += L(i, j) * L(i, j);
    }
    L(i, i) = std::sqrt(std::max(0.0, 1.0 - running_sum));
  }

  // The |dR/dL| part: empirically fitted against finite-difference
  // Jacobians (R^2 = 1 to numerical precision across P = 2..7) to be
  // sum_{i=1}^{P-1} (P - 1 - i) * log(L(i,i)) in 0-indexed terms.
  for(uword i = 1; i < P; i++) {
    double coef = (double)(P - 1 - i);
    if(coef > 0.0) {
      log_jac += coef * std::log(std::max(L(i, i), 1e-300));
    }
  }

  return log_jac;
};

//' @title Robbins-Monro adaptive proposal-window update
//' @description One diminishing-adaptation update of a scalar
//' Metropolis-Hastings proposal window; see the section header comment for
//' the references and the convergence argument. Operates in log-space so
//' the window stays strictly positive.
//' @param window Current (strictly positive) proposal window value.
//' @param acceptance_rate The realised acceptance rate this sweep, in
//' [0, 1]. May be a fraction over several components sharing one window
//' (e.g. the mean acceptance indicator across K clusters), not just a
//' single 0/1 draw - this is the standard "batched" extension of the
//' scalar Robbins-Monro process (Garthwaite et al., 2016, Section 3).
//' @param target_rate The target acceptance rate for this block: ~0.44 for
//' a scalar (1-D) random-walk update, ~0.234 for a block update that moves
//' several correlated dimensions at once (Roberts, Gelman & Gilks, 1997,
//' "Weak convergence and optimal scaling of random walk Metropolis
//' algorithms", Annals of Applied Probability 7(1)).
//' @param n The adaptation step index (e.g. the current MCMC iteration
//' within the burn-in window, 1-based). Larger n gives a smaller, more
//' conservative update, which is what makes the total adaptation finite.
//' @param step_scale Constant multiplying the 1/n^kappa step size.
//' @param kappa Decay exponent; must be in (0.5, 1] for the diminishing-
//' adaptation guarantee to apply. Default 0.6 follows Garthwaite et al.'s
//' recommendation.
//' @return The updated (still strictly positive) proposal window.
// [[Rcpp::export]]
double robbinsMonroUpdate(
  double window,
  double acceptance_rate,
  double target_rate,
  double n,
  double step_scale,
  double kappa
) {
  double step = step_scale / std::pow(std::max(n, 1.0), kappa);
  double log_window = std::log(window) + step * (acceptance_rate - target_rate);

  // Guard against a runaway update from a noisy early acceptance rate
  // (e.g. the very first sweep, where n = 1 gives the largest step size);
  // this keeps the window in a numerically sane range without otherwise
  // affecting the adaptation once it settles.
  log_window = std::min(std::max(log_window, -20.0), 20.0);

  return std::exp(log_window);
};
