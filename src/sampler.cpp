// # include <RcppArmadillo.h>
// # include <math.h>
// # include <string>
// # include <iostream>

# include "sampler.h"
# include "genericFunctions.h"
# include "pdfs.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;


// Parametrised class
sampler::sampler(
    arma::uword _K,
    arma::uword _B,
    arma::uvec _labels,
    arma::uvec _batch_vec,
    arma::vec _concentration,
    arma::mat _X,
    arma::uvec _fixed)
  {

    K = _K;
    B = _B;
    labels = _labels;
    batch_vec = _batch_vec;
    concentration = _concentration;
    X = _X;
    X_t = X.t();

    // Plausibly belongs in the MVN sampler. Used for selecting slices / columns
    // in the metropolis steps.
    KB_inds = linspace<uvec>(0, K - 1, K) * B;
    B_inds = linspace<uvec>(0, B - 1, B);

    // Dimensions
    N = X.n_rows;
    P = X.n_cols;

    // Missing-data bookkeeping (see sampler.h): X_raw/X_raw_t preserve
    // exactly what was passed in (including any NaN); items_to_augment
    // flags every row with at least one missing entry, generically, from
    // X_raw alone. mvnSamplerMixed additionally has censored/binary
    // entries to augment, so it recomputes (overwrites) this with its own
    // broader definition in its own constructor.
    X_raw = X;
    X_raw_t = X_t;
    std::vector<uword> rows_to_augment;
    for(uword n = 0; n < N; n++) {
      if(!X_t.col(n).is_finite()) {
        rows_to_augment.push_back(n);
      }
    }
    items_to_augment = conv_to<uvec>::from(rows_to_augment);

    // Class and batch populations
    N_k = zeros<uvec>(K);
    N_b = zeros<uvec>(B);

    // The batch labels won't ever change, so let's count them now
    for(uword b = 0; b < B; b++){
      N_b(b) = sum(batch_vec == b);
    }

    // Weights
    // double x, y;
    w = zeros<mat>(K, 1);

    // Log likelihood (individual and model)
    ll = zeros<vec>(K);
    likelihood = zeros<vec>(N);

    // Class members
    members.set_size(N, K);
    members.zeros();

    // The indices of the members of each batch in the dataset
    batch_ind.set_size(B);
    for(uword b = 0; b < B; b++) {
      batch_ind(b) = find(batch_vec == b);
    }

    // Fixed (semi-supervised) labels. A fixed vector of all zeroes gives the
    // fully unsupervised model.
    fixed = _fixed;
    N_fixed = arma::sum(fixed);
    uvec fixed_ind = find(fixed == 1);
    unfixed_ind = find(fixed == 0);

    // Allocation probability matrix. For fixed items this holds a one-hot
    // encoding of the known label; for unfixed items it is populated by
    // updateAllocation().
    alloc.set_size(N, K);
    alloc.zeros();
    for (auto& n : fixed_ind) {
      alloc(n, labels(n)) = 1.0;
    }

    // Interaction term and batch-specific weights default to entirely inert
    // (correctly-sized zero/identity structures) until/unless
    // initialiseInteraction()/initialiseBatchWeightPrior() are called from
    // a concrete sampler's constructor - this just guards against reading
    // uninitialised memory if a future sampler subclass forgets to opt in
    // explicitly.
    tau2_interaction = ones<vec>(P);
    gamma = zeros<cube>(P, K, B);
    gamma_count = zeros<uvec>(P);

    batch_coordinates = linspace<vec>(0, B - 1, B);
    gp_cov = eye<mat>(B, B);
    gp_chol = eye<mat>(B, B);
    w_batch = ones<mat>(B, K) / (double) K;
    eta_alr = zeros<mat>(B, (K > 0) ? K - 1 : 0);
    eta_count = zeros<uvec>((K > 0) ? K - 1 : 0);
    pp_mu = zeros<vec>((K > 0) ? K - 1 : 0);
    pp_tau2 = ones<vec>((K > 0) ? K - 1 : 0);
    gp_beta = zeros<vec>((K > 0) ? K - 1 : 0);
  };

// Functions required of all mixture models
// Function to update class / mixture weights
void sampler::updateWeights(){

  // Used locally as the posterior concentration
  double a = 0.0;

  for (uword k = 0; k < K; k++) {
    // Find how many labels have the value
    members.col(k) = labels == k;
    N_k(k) = sum(members.col(k));
  }

  if(weight_prior_type > 0) {
    // Batch-specific weights, either partial-pooling (exchangeable) or
    // GP-correlated - see the section header comment in sampler.h. w
    // itself is kept in sync as the N_b-weighted average across batches,
    // purely for any code that reports/reads the single global weight
    // vector, but updateAllocation() uses w_batch directly in this mode.
    if(weight_prior_type == 1) {
      updatePartialPoolingWeights();
    } else {
      updateGPWeights();
      if(sample_gp_hyperparameters) {
        gpHyperparameterMetropolis();
      }
    }
    return;
  }

  for (uword k = 0; k < K; k++) {
    // Update weights by sampling from a Gamma distribution
    a  = concentration(k) + (double) N_k(k);
    w(k) = randg( distr_param(a, 1.0) );
  }

  // Convert the cluster weights (previously gamma distributed) to Dirichlet
  // distributed by normalising (if K = 2 this is a Beta)
  w = w / accu(w);

};

// Sample the class allocations
void sampler::updateAllocation() {

  double u = 0.0, max_comp_prob = 0.0;
  uvec uniqueK;
  vec comp_prob(K);

  complete_likelihood = 0.0;
  observed_likelihood = 0.0;

  for(uword n = 0; n < N; n++){

    // The mixture-specific log likelihood for each observation in each class
    ll = itemLogLikelihood(X_t.col(n), batch_vec(n));

    // Update with weights: batch-specific under the GP-correlated weight
    // model, the single global vector otherwise (unchanged default path).
    if(weight_prior_type > 0) {
      comp_prob = ll + log(w_batch.row(batch_vec(n)).t());
    } else {
      comp_prob = ll + log(w);
    }

    // The observed-data log-likelihood contribution of item n. For a free
    // (unfixed) item, the label is a latent variable to be marginalised
    // over: log sum_k w_k f(x_n|theta_k) = log sum_k exp(comp_prob(k)),
    // computed here via the log-sum-exp trick for numerical stability, NOT
    // accu(comp_prob) (which would sum the log-densities across components
    // rather than marginalising over them). For a fixed (semi-supervised)
    // item the label is observed data, not something to marginalise over -
    // its contribution is the joint density at the known label,
    // w_{y_n} f(x_n|theta_{y_n}) = comp_prob(labels(n)), not a
    // marginalisation that also (wrongly) weighs in every other component
    // the item is known not to belong to.
    if(fixed(n) == 1) {
      observed_likelihood += comp_prob(labels(n));
    } else {
      max_comp_prob = max(comp_prob);
      observed_likelihood += max_comp_prob + log(accu(exp(comp_prob - max_comp_prob)));
    }

    // Handle overflow problems and then normalise to convert to probabilities
    comp_prob = exp(comp_prob - max_comp_prob);
    comp_prob = comp_prob / sum(comp_prob);

    // Prediction and update. The uniform draw is made unconditionally (even
    // for fixed items) so that the RNG draw sequence - and therefore the
    // sampled values for every other parameter - does not depend on which
    // items happen to be fixed.
    u = randu<double>( );

    if(fixed(n) == 0) {
      labels(n) = sum(u > cumsum(comp_prob));

      // The allocation probability for each class
      alloc.row(n) = comp_prob.t();
    }

    // Update the complete likelihood based on the new labelling
    complete_likelihood += ll(labels(n));
  }

  // Number of occupied components (used in BIC calculation)
  uniqueK = unique(labels);
  K_occ = uniqueK.n_elem;
};

// =============================================================================
// Batch x cluster interaction term (opt-in; see the header for the model).

void sampler::initialiseInteraction(
  bool _include_interaction,
  double _gamma_proposal_window,
  double _a_gamma,
  double _b_gamma
) {
  include_interaction = _include_interaction;
  gamma_proposal_window = _gamma_proposal_window;
  a_gamma = _a_gamma;
  b_gamma = _b_gamma;

  tau2_interaction = ones<vec>(P);
  gamma = zeros<cube>(P, K, B);
  gamma_count = zeros<uvec>(P);
};

// See the section header comment in sampler.h for the derivation and
// provenance of this identity.
arma::mat sampler::doubleCenterMatrix(arma::mat M) {
  vec row_mean = mean(M, 1);
  rowvec col_mean = mean(M, 0);
  double grand_mean = mean(row_mean);
  mat centered = M;
  centered.each_col() -= row_mean;
  centered.each_row() -= col_mean;
  centered += grand_mean;
  return centered;
};

// Draws gamma(p, ., .) from an (unconstrained) iid N(0, tau2_interaction(p))
// matrix, then projects it onto the sum-to-zero subspace via
// doubleCenterMatrix() - the projected draw is a proper (singular)
// zero-mean Gaussian supported on that subspace, and is the correct prior
// for a sum-to-zero-constrained interaction term (see sampler.h).
void sampler::sampleGammaPrior() {
  for(uword p = 0; p < P; p++) {
    mat raw(K, B);
    for(uword k = 0; k < K; k++) {
      for(uword b = 0; b < B; b++) {
        raw(k, b) = randn() * std::sqrt(tau2_interaction(p));
      }
    }
    mat centered = doubleCenterMatrix(raw);
    for(uword k = 0; k < K; k++) {
      for(uword b = 0; b < B; b++) {
        gamma(p, k, b) = centered(k, b);
      }
    }
  }
};

// Conjugate Gibbs update: gamma(p, ., .), confined to the sum-to-zero
// subspace, has exactly (K-1)(B-1) effective degrees of freedom per
// feature (not K*B - the remaining K*B - (K-1)(B-1) directions are always
// exactly zero by construction, and must not be counted as data
// contributing to tau2_interaction's posterior, or its evidence would be
// inflated and the posterior biased towards implausibly small values -
// checked directly: using K*B here measurably over-shrunk tau2_interaction
// relative to the (K-1)(B-1) version in a simulation with a genuine
// interaction effect). tau2_interaction(p) ~ InvGamma(a_gamma, b_gamma) is
// conjugate to this, giving posterior
// InvGamma(a_gamma + (K-1)(B-1)/2, b_gamma + 0.5 * sum gamma(p,.,.)^2).
void sampler::sampleTauInteractionPosterior() {
  double n_eff = (double) ((K - 1) * (B - 1));
  for(uword p = 0; p < P; p++) {
    double sum_sq = 0.0;
    for(uword k = 0; k < K; k++) {
      for(uword b = 0; b < B; b++) {
        sum_sq += gamma(p, k, b) * gamma(p, k, b);
      }
    }
    tau2_interaction(p) = rInvGamma(a_gamma + 0.5 * n_eff, b_gamma + 0.5 * sum_sq);
  }
};

// See the header comment: extra free parameters (beyond a concrete
// calcBIC()'s baseline "K_occ cluster blocks + B batch blocks" count) from
// the opt-in interaction term and/or batch-specific weight prior.
double sampler::structuralExtraBICParams() const {
  double extra = 0.0;

  if (include_interaction && K_occ > 0 && B > 0) {
    // gamma(p, ., .) has exactly (K_occ-1)(B-1) effective degrees of
    // freedom per feature once confined to the sum-to-zero subspace - see
    // sampleTauInteractionPosterior() above.
    extra += (double) (K_occ - 1) * (double) (B - 1) * (double) P;
  }

  if (weight_prior_type > 0 && K_occ > 0) {
    // Each concrete calcBIC() bakes in one weight parameter per occupied
    // component (K_occ, via its "1 +" cluster-block term), i.e. a single
    // shared weight vector. Under a batch-specific weight prior there is
    // instead a full (K_occ - 1)-free-coordinate simplex per batch, plus
    // the hierarchical prior's own hyperparameters.
    extra += (double) (K_occ - 1) * (double) B - (double) K_occ;

    if (weight_prior_type == 1) {
      // Partial pooling: population mean and variance for each of the
      // K_occ - 1 free ALR coordinates (pp_mu, pp_tau2).
      extra += 2.0 * (double) (K_occ - 1);
    } else if (weight_prior_type == 2) {
      // GP: one estimated intercept (gp_beta) per free ALR coordinate,
      // always - the GP analogue of partial pooling's per-coordinate
      // pp_mu (see the header for why this exists). tau2/length_scale are
      // each a single value SHARED across every coordinate (unlike
      // pp_mu/pp_tau2, which are per-coordinate), so they contribute only
      // 2 extra parameters in total, and only when actually estimated
      // rather than fixed by the user.
      extra += (double) (K_occ - 1);
      if (sample_gp_hyperparameters) {
        extra += 2.0;
      }
    }
  }

  return extra;
}

// =============================================================================
// Batch-specific mixing weights (opt-in; see the header for the model).

void sampler::initialiseBatchWeightPrior(
  arma::uword _weight_prior_type,
  arma::vec _batch_coordinates,
  double _gp_tau2,
  double _gp_length_scale,
  double _eta_proposal_window,
  bool _sample_gp_hyperparameters,
  double _gp_hyperparameter_proposal_window,
  double _pp_tau2_prior_shape,
  double _pp_tau2_prior_rate,
  double _pp_mu_prior_sd
) {
  weight_prior_type = _weight_prior_type;
  batch_coordinates = _batch_coordinates;
  gp_tau2 = _gp_tau2;
  gp_length_scale = _gp_length_scale;
  eta_proposal_window = _eta_proposal_window;
  sample_gp_hyperparameters = _sample_gp_hyperparameters;
  gp_hyperparameter_proposal_window = _gp_hyperparameter_proposal_window;
  pp_tau2_prior_shape = _pp_tau2_prior_shape;
  pp_tau2_prior_rate = _pp_tau2_prior_rate;
  pp_mu_prior_sd = _pp_mu_prior_sd;

  if (weight_prior_type == 2 && batch_coordinates.n_elem >= 2) {
    // Empirical-Bayes calibration of the length-scale prior to
    // batch_coordinates' own scale: an Inverse-Gamma "boundary-avoiding"
    // prior (Betancourt, 2020, "Robust Gaussian Process Modeling," Stan
    // case study; the same principle - discourage length scales shorter
    // than the finest resolvable spacing or longer than the whole domain -
    // underlies the penalised-complexity range prior of Fuglstad, Simpson,
    // Lindgren & Rue, 2019, JASA 114(525)) with its 1st/99th percentiles
    // placed AT those two bounds (invGammaQuantileMatch(), tail_prob =
    // 0.01), via a real quantile match rather than a log-normal
    // approximation - Inverse-Gamma's much faster-decaying right tail is
    // the actual reason it is the standard choice here (see the header):
    // a log-normal calibrated the same way was verified to still let the
    // sampler wander into, and get stuck at, length scales in the hundreds
    // to thousands on data where a single pair of batches happened to be
    // collected close together while the rest were spread out - Inverse-
    // Gamma does not have the tail mass left to permit that once its
    // quantiles are pinned at sensible bounds.
    //
    // The lower bound uses the MEDIAN gap between consecutive batches, not
    // the bare minimum: the minimum is a single order statistic and easily
    // dominated by one coincidentally-close pair (exactly what produced
    // the log-normal failure above), whereas the median reflects the
    // typical, "worth resolving" spacing the GP is actually meant to act
    // on. The upper bound is the full span of batch_coordinates, which has
    // no such single-pair fragility (Betancourt, 2020, uses the same
    // median-gap-to-full-range pairing).
    vec sorted_bc = arma::sort(batch_coordinates);
    vec gaps = arma::diff(sorted_bc);
    gaps = gaps.elem(arma::find(gaps > 1e-8)); // drop exact ties (zero gap)
    double lower_bound = gaps.n_elem > 0 ? arma::median(gaps) : 1.0;
    double range_bc = sorted_bc(sorted_bc.n_elem - 1) - sorted_bc(0);
    if (range_bc > 1e-8 && lower_bound > 1e-8 && range_bc > lower_bound) {
      invGammaQuantileMatch(lower_bound, range_bc, 0.01, gp_length_scale_prior_shape, gp_length_scale_prior_rate);
    }
    // If gp_length_scale itself was left at its own (possibly now
    // badly-scaled) default, re-centre it on the newly-calibrated prior's
    // MODE (rate / (shape + 1); more stable than the mean when shape is
    // close to 1, where the mean can itself be large or undefined) rather
    // than starting the chain somewhere the prior has just decided is
    // implausible.
    if (gp_length_scale == 1.0) {
      gp_length_scale = gp_length_scale_prior_rate / (gp_length_scale_prior_shape + 1.0);
    }
  }

  // tau2 (the marginal variance ON THE LOGIT SCALE) has no coordinate-unit
  // dependence to calibrate against, so it is not touched above - but note
  // the header default (gp_tau2_prior_shape = 2, gp_tau2_prior_rate = 4,
  // mean tau2 = 4) is itself already the weakly-informative logit-scale
  // variance Gelman, Jakulin, Pittau & Su (2008, Annals of Applied
  // Statistics 2(4), "A weakly informative default prior distribution for
  // logistic and other regression models") recommend as a generic default
  // for logistic-regression coefficients (SD around 2.5, i.e. variance
  // around 6.25) - widened from an earlier, tighter (mean 1) default that
  // predated this comparison.

  // No fixed proportional jitter is safe for every possible
  // batch_coordinates configuration: batches collected close together in
  // time elsewhere in the same dataset as a batch far from everyone else
  // can still make the Gram matrix arbitrarily ill-conditioned regardless
  // of tau2 - verified directly (condition number ~1.75e4, at this
  // package's own previously-fixed jitter proportion, for 10 realistically
  // irregularly-spaced batches whose closest pair happened to be far
  // closer together than their farthest) - see buildWellConditionedGPChol()
  // for the escalation this now uses instead.
  buildWellConditionedGPChol(gp_tau2, gp_length_scale, gp_cov, gp_chol, gp_jitter);

  // Start from equal weights in every batch (eta = 0 in every ALR
  // coordinate, and the partial-pooling population mean/variance at their
  // prior means); updateWeights() then updates this via MH/Gibbs from the
  // first sweep, exactly as the unconstrained w vector is only ever set
  // via the Gibbs step in updateWeights(), never a separate "sample from
  // prior" call.
  w_batch = ones<mat>(B, K) / (double) K;
  eta_alr = zeros<mat>(B, (K > 0) ? K - 1 : 0);
  eta_count = zeros<uvec>((K > 0) ? K - 1 : 0);
  pp_mu = zeros<vec>((K > 0) ? K - 1 : 0);
  pp_tau2 = ones<vec>((K > 0) ? K - 1 : 0) * (pp_tau2_prior_rate / std::max(pp_tau2_prior_shape - 1.0, 1e-6));
  // gp_beta plays the same role for "gp" that pp_mu plays for "partial
  // pooling" (see the header) - initialised at its own prior mean (0).
  gp_beta = zeros<vec>((K > 0) ? K - 1 : 0);
};

// Per-batch, per-cluster item counts under the current allocation, shared
// by every batch-weight-prior update.
arma::mat sampler::computeBatchClassCounts() {
  mat count_bk(B, K, fill::zeros);
  for(uword b = 0; b < B; b++) {
    for(uword k = 0; k < K; k++) {
      count_bk(b, k) = (double) accu((batch_vec == b) % (labels == k));
    }
  }
  return count_bk;
};

// Converts the current eta_alr (B x (K-1)) to per-batch simplex weights
// w_batch (B x K), and keeps the single global weight vector w in sync as
// the N_b-weighted average across batches, for any code that still reads
// w directly. Shared by every batch-weight-prior update.
void sampler::updateSimplexFromALR(arma::uword n_free) {
  vec N_b_d = conv_to<vec>::from(N_b);

  vec D = ones<vec>(B);
  for(uword j = 0; j < n_free; j++) {
    D += exp(eta_alr.col(j));
  }
  for(uword j = 0; j < n_free; j++) {
    w_batch.col(j) = exp(eta_alr.col(j)) / D;
  }
  w_batch.col(K - 1) = 1.0 / D;

  vec w_avg(K, fill::zeros);
  for(uword k = 0; k < K; k++) {
    w_avg(k) = accu(N_b_d % w_batch.col(k)) / (double) N;
  }
  w = w_avg;
};

// One sweep of the GP-structured batch-weight update: K - 1 independent
// ALR coordinates, each linked across the B batches by its own GP prior
// centred on its own estimated intercept gp_beta(j) (see the header for
// why a bare zero-mean GP is a needlessly restrictive special case), each
// updated by a block Metropolis-Hastings step (propose all B entries of
// that coordinate at once, accept/reject as a whole) holding every other
// coordinate fixed - the cyclic scheme documented alongside
// multinomialLogitGPLogKernel() (Ren, Du, Carin & Dunson, 2011; Linderman,
// Johnson & Adams, 2015). gp_beta(j) itself is then updated by an exact
// conjugate Gibbs step, generalised-least-squares rather than the simple
// sample mean partial pooling's analogous pp_mu update uses, since the
// deviations eta_{.,j} - gp_beta(j) are correlated (covariance gp_cov),
// not iid.
void sampler::updateGPWeights() {

  uword n_free = (K > 0) ? K - 1 : 0;
  mat count_bk = computeBatchClassCounts();
  vec N_b_d = conv_to<vec>::from(N_b);
  vec ones_B = ones<vec>(B);
  // Sigma^-1 * v == solve(L', solve(L, v)) exactly, via two triangular
  // solves against the Cholesky factor - never Sigma's explicit inverse;
  // see genericFunctions.h. Both solves below reuse this same L for every
  // free coordinate this sweep, since gp_chol only changes in
  // gpHyperparameterMetropolis().
  mat L = arma::trimatl(gp_chol);
  vec L_inv_ones = arma::solve(L, ones_B);

  vec eta_other_sum(B), eta_proposed(B), current_col(B);
  double proposed_score = 0.0, current_score = 0.0, u = 0.0, acceptance_prob = 0.0;

  for(uword j = 0; j < n_free; j++) {

    eta_other_sum.zeros();
    for(uword jp = 0; jp < n_free; jp++) {
      if(jp == j) continue;
      eta_other_sum += exp(eta_alr.col(jp));
    }

    current_col = eta_alr.col(j);
    eta_proposed = current_col + randn<vec>(B) * eta_proposal_window;

    proposed_score = multinomialLogitGPLogKernel(eta_proposed, eta_other_sum, count_bk.col(j), N_b_d, gp_chol, gp_beta(j));
    current_score = multinomialLogitGPLogKernel(current_col, eta_other_sum, count_bk.col(j), N_b_d, gp_chol, gp_beta(j));

    u = randu();
    acceptance_prob = std::min(1.0, std::exp(proposed_score - current_score));

    if(u < acceptance_prob) {
      eta_alr.col(j) = eta_proposed;
      eta_count(j)++;
    }

    // Gibbs update gp_beta(j) | eta_{.,j}, gp_cov (conjugate Normal-Normal
    // generalised least squares): the deviation model eta_{.,j} = beta_j *
    // 1 + w, w ~ N(0, gp_cov), combined with a N(0, pp_mu_prior_sd^2)
    // prior on beta_j, gives posterior precision 1' Sigma^-1 1 +
    // 1/pp_mu_prior_sd^2 and posterior mean (that precision)^-1 times
    // 1' Sigma^-1 eta_{.,j} - the ordinary GLS-with-a-prior formula
    // (Gelman et al., 2013, BDA3, Section 14.8), reducing to the simple
    // sample-mean update updatePartialPoolingWeights() uses for pp_mu
    // exactly when gp_cov is diagonal (no correlation to account for).
    // Both quadratic forms via the same triangular-solve identity as above
    // (1' Sigma^-1 1 = ||L^-1 1||^2, 1' Sigma^-1 eta = (L^-1 1)' (L^-1 eta)).
    vec L_inv_eta = arma::solve(L, eta_alr.col(j));
    double precision_beta = arma::dot(L_inv_ones, L_inv_ones) + 1.0 / (pp_mu_prior_sd * pp_mu_prior_sd);
    double post_var_beta = 1.0 / precision_beta;
    double post_mean_beta = post_var_beta * arma::dot(L_inv_ones, L_inv_eta);
    gp_beta(j) = post_mean_beta + randn() * std::sqrt(post_var_beta);
  }

  updateSimplexFromALR(n_free);
};

// One sweep of the partial-pooling (exchangeable) batch-weight update: for
// each free ALR coordinate j, eta_{b,j} | mu_j, tau2_j ~ N(mu_j, tau2_j)
// iid across batches b - no assumed order or distance between batches, in
// contrast to updateGPWeights(). eta_{.,j} itself is updated by a block
// Metropolis-Hastings step (as in the GP case, since the multinomial-logit
// likelihood isn't conjugate to anything), but with this simpler diagonal
// prior; mu_j and tau2_j are then updated by exact Gibbs steps, since
// (given eta_{.,j}) they form an ordinary conjugate Normal-Normal /
// Normal-InverseGamma hierarchical model (Gelman & Hill, 2007, ch. 12).
void sampler::updatePartialPoolingWeights() {

  uword n_free = (K > 0) ? K - 1 : 0;
  mat count_bk = computeBatchClassCounts();
  vec N_b_d = conv_to<vec>::from(N_b);

  vec eta_other_sum(B), eta_proposed(B), current_col(B);

  for(uword j = 0; j < n_free; j++) {

    eta_other_sum.zeros();
    for(uword jp = 0; jp < n_free; jp++) {
      if(jp == j) continue;
      eta_other_sum += exp(eta_alr.col(jp));
    }

    current_col = eta_alr.col(j);
    eta_proposed = current_col + randn<vec>(B) * eta_proposal_window;

    double log_lik_current = 0.0, log_lik_proposed = 0.0, D_b = 0.0;
    for(uword b = 0; b < B; b++) {
      D_b = 1.0 + std::exp(current_col(b)) + eta_other_sum(b);
      log_lik_current += count_bk(b, j) * current_col(b) - N_b_d(b) * std::log(D_b);
      D_b = 1.0 + std::exp(eta_proposed(b)) + eta_other_sum(b);
      log_lik_proposed += count_bk(b, j) * eta_proposed(b) - N_b_d(b) * std::log(D_b);
    }

    // Diagonal N(mu_j, tau2_j) prior - no matrix inversion needed, unlike
    // the GP case's full covariance.
    double prior_current = -0.5 * accu(square(current_col - pp_mu(j))) / pp_tau2(j);
    double prior_proposed = -0.5 * accu(square(eta_proposed - pp_mu(j))) / pp_tau2(j);

    double current_score = log_lik_current + prior_current;
    double proposed_score = log_lik_proposed + prior_proposed;

    double u = randu();
    double acceptance_prob = std::min(1.0, std::exp(proposed_score - current_score));

    if(u < acceptance_prob) {
      eta_alr.col(j) = eta_proposed;
      eta_count(j)++;
    }

    // Gibbs update mu_j | eta_{.,j}, tau2_j (conjugate Normal-Normal, prior
    // mu_j ~ N(0, pp_mu_prior_sd^2)).
    vec eta_col = eta_alr.col(j);
    double post_var = 1.0 / (1.0 / (pp_mu_prior_sd * pp_mu_prior_sd) + (double) B / pp_tau2(j));
    double post_mean = post_var * (accu(eta_col) / pp_tau2(j));
    pp_mu(j) = post_mean + randn() * std::sqrt(post_var);

    // Gibbs update tau2_j | eta_{.,j}, mu_j (conjugate InverseGamma).
    double sum_sq = accu(square(eta_col - pp_mu(j)));
    pp_tau2(j) = rInvGamma(pp_tau2_prior_shape + 0.5 * (double) B, pp_tau2_prior_rate + 0.5 * sum_sq);
  }

  updateSimplexFromALR(n_free);
};

// Builds a Matern-3/2 GP covariance matrix and its LOWER-triangular
// Cholesky factor for a candidate (tau2, length_scale) pair, escalating
// the diagonal jitter geometrically from a tau2-scaled starting point
// until chol() succeeds - the modern standard safeguard for a covariance
// matrix whose conditioning depends on the (here, potentially very
// irregular) input locations, not just on tau2/length_scale in isolation
// (GPyTorch's/GPflow's escalating default jitter; Gardner, Pleiss, Bindel,
// Weinberger & Wilson, 2018, "GPyTorch," NeurIPS). A single fixed
// proportion of tau2 (this package's first attempt at this fix) is not
// safe for every batch_coordinates configuration: verified directly to
// still leave a condition number of ~1.75e4 for a realistic set of 10
// irregularly-spaced batches whose closest pair is far closer together
// than its farthest. Deliberately factorises via chol(), not inv_sympd():
// every consumer (multinomialLogitGPLogKernel(), the whitening step in
// gpHyperparameterMetropolis() below) only ever needs a triangular solve
// against this factor, never the covariance's explicit inverse - see
// genericFunctions.h for why that is the numerically-preferred approach.
void sampler::buildWellConditionedGPChol(
  double tau2, double length_scale,
  arma::mat& cov, arma::mat& chol_factor, double& jitter_used
) {
  double jitter = std::max(1e-6, tau2 * 1e-4);
  const uword max_attempts = 12; // 1e-4 * 10^12 saturates long before this

  for (uword attempt = 0; attempt < max_attempts; attempt++) {
    cov = maternKernel32(batch_coordinates, tau2, length_scale, jitter);
    if (arma::chol(chol_factor, cov, "lower")) {
      jitter_used = jitter;
      return;
    }
    jitter *= 10.0;
  }

  // Every escalation attempt failed (only possible with a pathological
  // batch_coordinates input, e.g. hundreds of exactly-duplicated
  // timestamps) - fall back to the final, most-regularised attempt's
  // covariance and accept whatever chol() manages rather than leaving
  // chol_factor stale/uninitialised; a covariance this heavily jittered is
  // already so close to jitter * I that chol() succeeding is not actually
  // in doubt in practice.
  jitter_used = jitter / 10.0;
  cov = maternKernel32(batch_coordinates, tau2, length_scale, jitter_used);
  arma::chol(chol_factor, cov, "lower");
};

// Metropolis-Hastings update for the GP marginal variance (tau2) and
// length scale, proposed jointly as a symmetric random walk on their logs
// (keeps both strictly positive), together with a matching update to
// eta_alr via a NON-CENTRED/WHITENED reparameterisation used for this step
// only (the rest of the sampler stays in the ordinary "centred" eta_alr
// representation - see updateGPWeights()).
//
// Why: proposing new hyperparameters while eta_alr stays rigidly fixed
// (this function's previous approach, and the textbook-obvious one) directly
// couples the acceptance ratio to how "surprised" the CURRENT eta_alr is
// under the new covariance - the classic Neal's-funnel pathology of jointly
// modelling a hierarchical variance/length-scale parameter and the latent
// values it governs in their natural, entangled form (Neal, 2003, "Slice
// sampling," Annals of Statistics, Section 8.1; Papaspiliopoulos, Roberts &
// Skold, 2007, "A general framework for the parametrization of hierarchical
// models," Statistical Science; Betancourt & Girolami, 2015, "Hamiltonian
// Monte Carlo for Hierarchical Models"). Betancourt (2020, "Robust Gaussian
// Process Modeling," Stan case study) demonstrates the same pathology for a
// GP's length-scale specifically, and its standard fix: work instead with
// the ancillary, hyperparameter-FREE quantity z = L_current^-1 (eta - beta)
// (a triangular solve, not an inversion), propose new hyperparameters, then
// re-map the SAME z through the PROPOSED covariance's Cholesky factor to get
// the eta a well-behaved joint move implies. z's own N(0, I) prior density
// is identical before and after (z itself never changed, only what it maps
// to), so it cancels exactly and never needs to be evaluated at all - the
// acceptance ratio is left with only the (data) likelihood difference
// between the current and implied eta, plus the ordinary priors on
// tau2/length_scale themselves. This is also why gpHyperparameterLogKernel()
// - which evaluated a GP prior density under two different covariance
// matrices directly, needing their log-determinants - no longer exists:
// centring the comparison on the data likelihood instead removes the need
// for it, and for gp_cov_inv, entirely.
void sampler::gpHyperparameterMetropolis() {

  uword n_free = (K > 0) ? K - 1 : 0;
  mat count_bk = computeBatchClassCounts();
  vec N_b_d = conv_to<vec>::from(N_b);

  double log_tau2_current = std::log(gp_tau2),
    log_length_scale_current = std::log(gp_length_scale);

  double log_tau2_proposed = log_tau2_current + randn() * gp_hyperparameter_proposal_window,
    log_length_scale_proposed = log_length_scale_current + randn() * gp_hyperparameter_proposal_window;

  double tau2_proposed = std::exp(log_tau2_proposed),
    length_scale_proposed = std::exp(log_length_scale_proposed);

  // buildWellConditionedGPChol() re-derives whatever jitter this particular
  // (tau2, length_scale) pair actually needs (see its own documentation) -
  // simply tracking the previous jitter, as an earlier version of this
  // function did, is not sufficient on its own, since escalation is
  // driven by the *matrix's* conditioning, which depends on length_scale
  // and batch_coordinates too, not tau2 alone.
  mat gp_cov_proposed, gp_chol_proposed;
  double jitter_proposed = 0.0;
  buildWellConditionedGPChol(tau2_proposed, length_scale_proposed, gp_cov_proposed, gp_chol_proposed, jitter_proposed);

  mat L_current = arma::trimatl(gp_chol);
  mat L_proposed = arma::trimatl(gp_chol_proposed);
  mat eta_implied(B, n_free);
  for (uword j = 0; j < n_free; j++) {
    vec z_j = arma::solve(L_current, eta_alr.col(j) - gp_beta(j));
    eta_implied.col(j) = gp_beta(j) + L_proposed * z_j;
  }

  double log_lik_current = 0.0, log_lik_implied = 0.0;
  for (uword j = 0; j < n_free; j++) {
    vec eta_other_sum_cur(B, fill::zeros), eta_other_sum_imp(B, fill::zeros);
    for (uword jp = 0; jp < n_free; jp++) {
      if (jp == j) continue;
      eta_other_sum_cur += exp(eta_alr.col(jp));
      eta_other_sum_imp += exp(eta_implied.col(jp));
    }
    log_lik_current += multinomialLogitLogLik(eta_alr.col(j), eta_other_sum_cur, count_bk.col(j), N_b_d);
    log_lik_implied += multinomialLogitLogLik(eta_implied.col(j), eta_other_sum_imp, count_bk.col(j), N_b_d);
  }

  // Both tau2 and length_scale now have Inverse-Gamma priors expressed
  // directly in their own (original, not log-) scale, so sampling either
  // via a symmetric walk on its log needs the usual explicit
  // change-of-variables Jacobian (+log(proposed value)) in the acceptance
  // ratio - length_scale's prior switched from a log-normal (which built
  // the Jacobian into the density term itself, needing none added
  // separately - see the earlier version of this function) to Inverse-
  // Gamma specifically for its far-faster-decaying right tail; see the
  // header for why that tail behaviour is what actually matters here.
  double proposed_score = log_lik_implied
    + invGammaLogLikelihood(tau2_proposed, gp_tau2_prior_shape, gp_tau2_prior_rate)
    + invGammaLogLikelihood(length_scale_proposed, gp_length_scale_prior_shape, gp_length_scale_prior_rate)
    + log_tau2_proposed + log_length_scale_proposed;

  double current_score = log_lik_current
    + invGammaLogLikelihood(gp_tau2, gp_tau2_prior_shape, gp_tau2_prior_rate)
    + invGammaLogLikelihood(gp_length_scale, gp_length_scale_prior_shape, gp_length_scale_prior_rate)
    + log_tau2_current + log_length_scale_current;

  double u = randu();
  double acceptance_prob = std::min(1.0, std::exp(proposed_score - current_score));

  if(u < acceptance_prob) {
    gp_tau2 = tau2_proposed;
    gp_length_scale = length_scale_proposed;
    gp_jitter = jitter_proposed;
    gp_cov = gp_cov_proposed;
    gp_chol = gp_chol_proposed;
    eta_alr = eta_implied; // eta moves together with the hyperparameters - the whole point
    updateSimplexFromALR(n_free); // w_batch must reflect the now-changed eta_alr
    gp_hyperparameter_count++;
  }
};
