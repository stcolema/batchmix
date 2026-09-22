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
    gp_cov_inv = eye<mat>(B, B);
    w_batch = ones<mat>(B, K) / (double) K;
    eta_alr = zeros<mat>(B, (K > 0) ? K - 1 : 0);
    eta_count = zeros<uvec>((K > 0) ? K - 1 : 0);
    pp_mu = zeros<vec>((K > 0) ? K - 1 : 0);
    pp_tau2 = ones<vec>((K > 0) ? K - 1 : 0);
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

    // The observed (marginal, label-free) log-likelihood of item n is
    // log sum_k exp(comp_prob(k)), computed here via the log-sum-exp trick
    // for numerical stability, NOT accu(comp_prob) (which would sum the
    // log-densities across components rather than marginalising over them).
    max_comp_prob = max(comp_prob);
    observed_likelihood += max_comp_prob + log(accu(exp(comp_prob - max_comp_prob)));

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

  // A fixed, tiny absolute jitter (e.g. 1e-6) is nowhere near enough once
  // length_scale is not small relative to the spacing of batch_coordinates
  // (exactly the "batches are smoothly correlated over time" regime this
  // feature exists for): neighbouring batches then have near-identical
  // rows/columns, the smallest eigenvalue of gp_cov collapses towards the
  // jitter floor, and gp_cov_inv's corresponding eigenvalue explodes -
  // verified directly (condition number ~1.7e7 with tau2 = 5, length_scale
  // = 3 over 8 unit-spaced batches, jitter = 1e-6), which silently corrupts
  // the quadratic form in every *LogKernel/gpHyperparameterLogKernel call
  // using gp_cov_inv, in one observed case flipping which of two eta
  // vectors the GP prior term favours. Scaling the jitter with tau2 (the
  // matrix's own diagonal scale) keeps the condition number bounded
  // regardless of length_scale; 1e-4 * tau2 was verified to bring the
  // condition number for the case above down to ~2.2e3 while leaving the
  // resulting log-kernel differences essentially unchanged from their
  // well-conditioned limit. Only meaningful for weight_prior_type == 2
  // ("gp"), but harmless to always compute.
  gp_jitter = std::max(gp_jitter, gp_tau2 * 1e-4);
  gp_cov = squaredExponentialKernel(batch_coordinates, gp_tau2, gp_length_scale, gp_jitter);
  gp_cov_inv = inv_sympd(gp_cov);

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
// ALR coordinates, each linked across the B batches by its own GP prior,
// each updated by a block Metropolis-Hastings step (propose all B entries
// of that coordinate at once, accept/reject as a whole) holding every
// other coordinate fixed - the cyclic scheme documented alongside
// multinomialLogitGPLogKernel() (Ren, Du, Carin & Dunson, 2011; Linderman,
// Johnson & Adams, 2015).
void sampler::updateGPWeights() {

  uword n_free = (K > 0) ? K - 1 : 0;
  mat count_bk = computeBatchClassCounts();
  vec N_b_d = conv_to<vec>::from(N_b);

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

    proposed_score = multinomialLogitGPLogKernel(eta_proposed, eta_other_sum, count_bk.col(j), N_b_d, gp_cov, gp_cov_inv);
    current_score = multinomialLogitGPLogKernel(current_col, eta_other_sum, count_bk.col(j), N_b_d, gp_cov, gp_cov_inv);

    u = randu();
    acceptance_prob = std::min(1.0, std::exp(proposed_score - current_score));

    if(u < acceptance_prob) {
      eta_alr.col(j) = eta_proposed;
      eta_count(j)++;
    }
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

// The (unnormalised, additive-constant-free) log-density of a zero-mean
// P-variate... here B-variate... Gaussian with covariance `cov`, evaluated
// jointly at every free ALR coordinate's current value. Unlike the
// eta-coordinate MH step above (where gp_cov is held fixed and its log
// |.| term cancels in the ratio), a step that changes the GP
// hyperparameters themselves changes gp_cov, so that determinant term does
// NOT cancel and must be included here.
double sampler::gpHyperparameterLogKernel(arma::mat cov, arma::mat cov_inv) {
  uword n_free = (K > 0) ? K - 1 : 0;
  double log_det_val = log_det(cov).real();
  double score = 0.0;
  for(uword j = 0; j < n_free; j++) {
    score += -0.5 * (log_det_val + as_scalar(eta_alr.col(j).t() * cov_inv * eta_alr.col(j)));
  }
  return score;
};

// Metropolis-Hastings update for the GP marginal variance (tau2) and
// length scale, proposed jointly as a symmetric random walk on their logs
// (keeps both strictly positive). Proposing on the log scale while the
// prior is expressed in the original scale requires the usual
// change-of-variables Jacobian correction (+ log(proposed) - log(current)
// for each log-transformed parameter) in the acceptance ratio.
void sampler::gpHyperparameterMetropolis() {

  double log_tau2_current = std::log(gp_tau2),
    log_length_scale_current = std::log(gp_length_scale);

  double log_tau2_proposed = log_tau2_current + randn() * gp_hyperparameter_proposal_window,
    log_length_scale_proposed = log_length_scale_current + randn() * gp_hyperparameter_proposal_window;

  double tau2_proposed = std::exp(log_tau2_proposed),
    length_scale_proposed = std::exp(log_length_scale_proposed);

  // Jitter must track the proposed tau2 (see initialiseBatchWeightPrior()
  // for why a fixed jitter is unsafe): otherwise, once tau2 grows during
  // sampling, a jitter sized for the old (smaller) tau2 can again be too
  // small relative to the new diagonal scale.
  double jitter_proposed = std::max(gp_jitter, tau2_proposed * 1e-4);
  mat gp_cov_proposed = squaredExponentialKernel(batch_coordinates, tau2_proposed, length_scale_proposed, jitter_proposed);
  mat gp_cov_inv_proposed = inv_sympd(gp_cov_proposed);

  // tau2's prior (InvGamma) is expressed directly in tau2-space, so
  // sampling it via a symmetric walk on log(tau2) needs an explicit
  // +log(tau2) Jacobian term. length_scale's prior below is instead
  // written directly as a Normal kernel on log(length_scale) - i.e. a
  // LogNormal(mean, sd^2) prior on length_scale itself - which already
  // *is* "LogNormal density + Jacobian" combined (the 1/length_scale
  // factor in the LogNormal density and the Jacobian factor of
  // length_scale exactly cancel), so no separate Jacobian term is added
  // for it here.
  double proposed_score = gpHyperparameterLogKernel(gp_cov_proposed, gp_cov_inv_proposed)
    + invGammaLogLikelihood(tau2_proposed, gp_tau2_prior_shape, gp_tau2_prior_rate)
    - std::pow(log_length_scale_proposed - gp_length_scale_prior_mean, 2.0) / (2.0 * gp_length_scale_prior_sd * gp_length_scale_prior_sd)
    + log_tau2_proposed;

  double current_score = gpHyperparameterLogKernel(gp_cov, gp_cov_inv)
    + invGammaLogLikelihood(gp_tau2, gp_tau2_prior_shape, gp_tau2_prior_rate)
    - std::pow(log_length_scale_current - gp_length_scale_prior_mean, 2.0) / (2.0 * gp_length_scale_prior_sd * gp_length_scale_prior_sd)
    + log_tau2_current;

  double u = randu();
  double acceptance_prob = std::min(1.0, std::exp(proposed_score - current_score));

  if(u < acceptance_prob) {
    gp_tau2 = tau2_proposed;
    gp_length_scale = length_scale_proposed;
    gp_jitter = jitter_proposed;
    gp_cov = gp_cov_proposed;
    gp_cov_inv = gp_cov_inv_proposed;
    gp_hyperparameter_count++;
  }
};
