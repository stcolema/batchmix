// generic_functions.h
// =============================================================================
// include guard
#ifndef GENFUN_H
#define GENFUN_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>

using namespace arma ;

// =============================================================================
// a set of functions used in a few different places

//' @title Propose new non-negative value
//' @description Propose new non-negative for sampling.
//' @param x Current value to be proposed
//' @param window The proposal window
//' @return new double
double proposeNewNonNegativeValue(double x, double window, 
                                  bool use_log_norm = false,
                                  double tolerance = 1e-12
);

//' @title The Inverse Gamma Distribution
//' @description Random generation from the inverse Gamma distribution.
//' @param shape Shape parameter.
//' @param rate Rate parameter.
//' @return Sample from invGamma(shape, rate).
double rInvGamma(double shape, double rate);

//' @title The Inverse Gamma Distribution
//' @description Random generation from the inverse Gamma distribution.
//' @param N Number of samples to draw.
//' @param shape Shape parameter.
//' @param rate Rate parameter.
//' @return Sample from invGamma(shape, rate).
arma::vec rInvGamma(uword N, double shape, double rate);

//' @title The Gamma Distribution
//' @description Random generation from the Gamma distribution.
//' @param shape Shape parameter.
//' @param rate Rate parameter.
//' @return Sample from Gamma(shape, rate).
double rGamma(double shape, double rate);

//' @title The Gamma Distribution
//' @description Random generation from the Gamma distribution.
//' @param N Number of samples to draw.
//' @param shape Shape parameter.
//' @param rate Rate parameter.
//' @return N samples from Gamma(shape, rate).
arma::vec rGamma(uword N, double shape, double rate);

//' @title The Beta Distribution
//' @description Random generation from the Beta distribution.
//' See https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions.
//' Samples from a Beta distribution based using two independent gamma
//' distributions.
//' @param a Shape parameter.
//' @param b Shape parameter.
//' @return Sample from Beta(a, b).
double rBeta(double a, double b);

//' @title The Beta Distribution
//' @description Random generation from the Beta distribution.
//' See https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions.
//' Samples from a Beta distribution based using two independent gamma
//' distributions.
//' @param n The number of samples to draw.
//' @param a Shape parameter.
//' @param b Shape parameter.
//' @return Sample from Beta(a, b).
arma::vec rBeta(arma::uword n, double a, double b);

//' @title The Log-Normal Distribution
//' @description Random generation from the log-Normal distribution.
//' @param mu mean parameter.
//' @param sd standard deviation parameter.
//' @return Sample from log-Normal(mu, sd^2).
double rLogNormal(double mu, double sd);

//' @title The Log-Normal Distribution
//' @description Random generation from the log-Normal distribution.
//' @param N positive integer - the number of samples drawn.
//' @param mu mean parameter.
//' @param sd standard deviation parameter.
//' @return N samples from log-Normal(mu, sd^2)
arma::vec rLogNormal(arma::uword N, double mu, double sd);

//' @title Metropolis acceptance step
//' @description Given a probaility, randomly accepts by sampling from a uniform 
//' distribution.
//' @param acceptance_prob Double between 0 and 1.
//' @return Boolean indicating acceptance.
bool metropolisAcceptanceStep(double acceptance_prob);


//' title Accept proposal
//' description Determines if a proposal is accepted given a log ratio of scores 
//' for the proposed and original values.
//' param proposed_model_score Score in the posterior kernel for the proposed 
//' parameter value
//' param current_model_score Score in the posterior kernel for the current 
//' parameter value
//' return Boolean indicating acceptance.
bool acceptProposal(double proposed_model_score, double current_model_score);

//' @title Sample mean
//' @description calculate the sample mean of a matrix X.
//' @param X Matrix
//' @return Vector of the column means of X.
vec sampleMean(arma::mat X);

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
);

//' @title The LKJ distribution
//' @description Random generation of a P x P correlation matrix from the
//' LKJ(eta) distribution (Lewandowski, Kurowicka & Joe, 2009). Implemented
//' by rejection sampling: the free upper-triangular entries are drawn
//' i.i.d. Uniform(-1, 1) and the candidate is rejected if it is not
//' positive definite (this alone samples LKJ(eta = 1), i.e. the uniform
//' distribution over the space of correlation matrices); a further
//' accept/reject step with probability det(R)^(eta - 1) re-weights this to
//' general eta, which is a valid rejection sampler whenever eta >= 1
//' because det(R) <= 1 for any correlation matrix (Hadamard's inequality),
//' so det(R)^(eta - 1) <= 1 is a proper acceptance probability. For
//' eta < 1 the LKJ density is unbounded as R approaches singularity, so
//' this construction does not apply and is not supported here.
//' @param P Dimension of the correlation matrix.
//' @param eta Concentration parameter, eta >= 1. eta = 1 is uniform over
//' correlation matrices; eta > 1 concentrates mass towards the identity
//' (i.e. shrinks correlations towards 0).
//' @return A P x P correlation matrix sampled from LKJ(eta).
arma::mat sampleLKJCorrelationMatrix(arma::uword P, double eta);

//' @title The truncated Normal distribution
//' @description Random generation from a univariate Normal(mean, sd^2)
//' distribution truncated to the interval (lower, upper), via the inverse
//' CDF method. This is numerically reliable for truncation points within a
//' few standard deviations of the mean, which covers the overwhelming
//' majority of practical missing/censored-data and probit-augmentation use
//' cases; it is not the numerically robust choice for truncation deep in
//' the tail (e.g. beyond roughly 8-10 SDs), where the CDF/inverse-CDF
//' values saturate to 0/1 in floating point. For that regime, use an
//' exponential-tilting rejection sampler instead (Robert, 1995,
//' "Simulation of truncated normal variables", Statistics and Computing
//' 5(2)) - see also the 'truncnorm' or 'TruncatedNormal' R packages, which
//' implement this.
//' @param mean Mean of the untruncated Normal distribution.
//' @param sd Standard deviation of the untruncated Normal distribution.
//' @param lower Lower truncation bound (use -arma::datum::inf for none).
//' @param upper Upper truncation bound (use arma::datum::inf for none).
//' @return A draw from Normal(mean, sd^2) truncated to (lower, upper).
double rTruncNorm(double mean, double sd, double lower, double upper);

//' @title Squared-exponential covariance kernel
//' @description Builds the covariance matrix of a zero-mean Gaussian
//' process with a squared-exponential (RBF) kernel evaluated at a set of
//' 1-D locations, k(x, x') = tau2 * exp(-(x - x')^2 / (2 * length_scale^2)),
//' plus a small jitter added to the diagonal for numerical positive
//' definiteness. Intended use here: x is the (ordered, e.g. time-indexed)
//' batch index, so that this encodes "batch i's parameter is correlated
//' with batch i - 1 and i + 1, decaying with |i - j| / length_scale" -
//' see the design note in mvnSamplerMixed.h / the package documentation
//' for the batch-dependent-weights model this supports.
//' @param x Vector of 1-D locations (e.g. batch indices or timestamps).
//' @param tau2 Marginal variance (amount of variation away from the mean
//' function).
//' @param length_scale How many units of x the correlation persists over;
//' larger values give smoother (more strongly correlated) functions.
//' @param jitter Small value added to the diagonal for numerical
//' stability (e.g. 1e-6).
//' @return The covariance matrix, length(x) x length(x).
arma::mat squaredExponentialKernel(arma::vec x, double tau2, double length_scale, double jitter);

//' @title Multinomial-logit Gaussian process log-kernel
//' @description The unnormalised log-posterior-kernel for one additive
//' log-ratio (ALR) coordinate of a set of batch-dependent multinomial
//' class-weight vectors, under a Gaussian process prior over the ordered
//' batch index (Aitchison, 1982, for the ALR transform; Ren, Du, Carin &
//' Dunson, 2011, "The Logistic Stick-Breaking Process," JMLR, and
//' Linderman, Johnson & Adams, 2015, "Dependent Multinomial Models Made
//' Easy," NeurIPS, for GP/AR-linked multinomial weights - this is the
//' Metropolis-Hastings-compatible route described there, rather than
//' their more efficient but considerably more involved Polya-Gamma Gibbs
//' sampler). Combines the GP log-prior for this ALR coordinate across all
//' batches with the multinomial-logit log-likelihood contribution that
//' coordinate makes to every batch's observed class counts.
//' @param eta B-vector, the j-th ALR coordinate (log(w_{b,j} / w_{b,K}))
//' for each batch b = 1, ..., B.
//' @param eta_other_sum B-vector, sum_{j' != j} exp(eta_{b,j'}) for each
//' batch b (i.e. everything other than coordinate j and the pivot
//' category K that enters the softmax denominator); pass a vector of
//' zeros if K = 2 (only one free coordinate, no other categories).
//' @param class_counts_j B-vector, the number of items in each batch
//' currently allocated to class j (the numerator category for this ALR
//' coordinate).
//' @param class_counts_total B-vector, the total number of items in each
//' batch.
//' @param gp_cov The B x B Gaussian process covariance matrix for this
//' coordinate (see squaredExponentialKernel()).
//' @param gp_cov_inv The (precomputed) inverse of gp_cov.
//' @return The unnormalised log-posterior-kernel value for eta.
double multinomialLogitGPLogKernel(
  arma::vec eta,
  arma::vec eta_other_sum,
  arma::vec class_counts_j,
  arma::vec class_counts_total,
  arma::mat gp_cov,
  arma::mat gp_cov_inv
);

// =============================================================================
// Cholesky-factor reparameterisation of a correlation matrix, via canonical
// partial correlations, for random-walk Metropolis proposals on R.
//
// Random-walk proposals directly on R's pairwise entries (perturb, then
// reject if not positive definite) do not respect the geometry of the
// space of correlation matrices, and mix increasingly poorly as P grows
// and as more entries are simultaneously correlated (Stan Reference
// Manual, "Cholesky Factors of Correlation Matrices"; this is exactly why
// Stan's own lkj_corr_cholesky works in this parameterisation rather than
// on R directly). The fix used here: parameterise R via its Cholesky
// factor L (R = LL'), whose rows are unit vectors; each row i has i - 1
// free entries Z(i, j) in (-1, 1) (0-indexed; row 0 is trivially (1)),
// built up via a sequential "remaining L2-norm budget" construction
// (Joe, 2006, "Generating random correlation matrices based on partial
// correlations", J. Multivariate Analysis - the same construction
// underlying the LKJ paper's own "onion method"). Unlike raw pairwise
// correlations, ANY Z(i,j) in (-1,1) - independently, with no joint
// constraint - yields a valid positive definite R: no rejection is ever
// needed. See buildCorrelationCholeskyFromZ() for the construction and
// logJacobianZToR() for its (independently, numerically verified -
// against finite-difference Jacobians, to within 1e-9, at P = 2..7 -
// rather than merely recalled) Jacobian.

//' @title Build a correlation-matrix Cholesky factor from partial
//' correlations
//' @description Constructs the lower-triangular Cholesky factor L (unit-
//' norm rows, so that R = LL' is always a valid correlation matrix) from
//' free parameters Z(i, j), 0 <= j < i < P, each in (-1, 1). See the
//' section header comment for the construction and its provenance; any Z
//' in this domain yields a valid PD correlation matrix, so this never
//' needs to reject.
//' @param Z A P x P matrix; only the strictly-lower-triangular entries
//' Z(i,j), j < i, are used.
//' @param P The dimension.
//' @return The P x P lower-triangular Cholesky factor L.
arma::mat buildCorrelationCholeskyFromZ(arma::mat Z, arma::uword P);

//' @title Invert buildCorrelationCholeskyFromZ()
//' @description Recovers the partial correlations Z from a valid
//' correlation-matrix Cholesky factor L (lower-triangular, unit-norm
//' rows).
//' @param L The P x P Cholesky factor.
//' @param P The dimension.
//' @return The P x P matrix Z (only entries below the diagonal are
//' meaningful).
arma::mat choleskyToPartialCorrelations(arma::mat L, arma::uword P);

//' @title Jacobian of the partial-correlation to correlation-matrix map
//' @description log|dR/dZ|, the log-Jacobian determinant of the map from
//' the free partial correlations Z (see buildCorrelationCholeskyFromZ())
//' to the free entries of R = LL'. Needed to correctly do a Metropolis
//' random walk in an unconstrained transform of Z while targeting the
//' correct density for R. Verified against finite-difference Jacobians;
//' see the section header comment.
//' @param Z A P x P matrix; only the strictly-lower-triangular entries
//' are used (as in buildCorrelationCholeskyFromZ()).
//' @param P The dimension.
//' @return The log-Jacobian determinant, a scalar.
double logJacobianZToR(arma::mat Z, arma::uword P);

// =============================================================================
// Auto-tuning: Robbins-Monro adaptive proposal windows.
//
// Every Metropolis-Hastings proposal window in the package (mu, m, S,
// cov/r/sigma, t_df, and the GP-weight/interaction windows added alongside
// this) was previously a fixed constant the user had to hand-tune by trial
// and error. This implements diminishing-adaptation MH tuning (Roberts &
// Rosenthal, 2009, "Examples of Adaptive MCMC", Journal of Computational and
// Graphical Statistics 18(2); Garthwaite, Fan & Sisson, 2016, "Adaptive
// optimal scaling of Metropolis-Hastings algorithms using the Robbins-Monro
// process", Communications in Statistics - Theory and Methods 45(17)): the
// log-window is nudged towards whatever value would have produced a target
// acceptance rate, with a step size that shrinks as 1/n^kappa. The shrinking
// step size is what makes this valid - it guarantees the total amount of
// adaptation is finite, so the chain still converges to the correct
// stationary distribution once adaptation is switched off (this package
// freezes adaptation at the end of a user-supplied burn-in window, rather
// than continuing to adapt forever, as the simplest sufficient condition).

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
double robbinsMonroUpdate(
  double window,
  double acceptance_rate,
  double target_rate,
  double n,
  double step_scale = 1.0,
  double kappa = 0.6
);

#endif /* GENFUN_H */
