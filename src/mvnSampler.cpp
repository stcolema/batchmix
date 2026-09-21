// mvnSampler.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "pdfs.h"
# include "sampler.h"
# include "mvnSampler.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// mvnSampler class

mvnSampler::mvnSampler(                           
  arma::uword _K,
  arma::uword _B,
  double _mu_proposal_window,
  double _cov_proposal_window,
  double _m_proposal_window,
  double _S_proposal_window,
  arma::uvec _labels,
  arma::uvec _batch_vec,
  arma::vec _concentration,
  arma::mat _X,
  double _m_scale,
  double _rho,
  double _theta,
  bool _sample_m_scale
) : sampler(_K,
_B,
_labels,
_batch_vec,
_concentration,
_X) {
  
  rowvec X_min = min(X), X_max = max(X);
  mat global_cov = arma::cov(X);
  
  n_param_cluster = 1 + P + P * (P + 1) * 0.5;
  n_param_batch = 2 * P;
  nu = P + 2;
  
  // Mean
  mat mean_mat = mean(_X, 0).t();
  xi = mean_mat.col(0);
  
  // Empirical Bayes fora diagonal covariance matrix
  mat scale_param = _X.each_row() - xi.t();
  vec diag_entries(P);
  // double scale_entry = accu(scale_param % scale_param, 0) / (N * std::pow(K, 1.0 / (double) P));
  
  double scale_entry = (accu(global_cov.diag()) / P) / std::pow(K, 2.0 / (double) P);
  
  diag_entries.fill(scale_entry);
  scale = diagmat( diag_entries );
  
  scale = global_cov / std::pow(K, 2.0 / (double) P);
  
  // The mean of the prior distribution for the batch shift, m, parameter
  delta_2 = 0.0;
  lambda_2 = _m_scale;
  m_scale = _m_scale;
  
  // Prior precision is the inverse of something on the scale of 1/10 the global 
  // covariance
  delta_2 = accu(global_cov.diag()) / (double) P ;
  batch_shift_prior_precision = 1.0 / (delta_2 * lambda_2);
  
  t = 1.0 / ((accu(global_cov.diag()) / P ) * m_scale);
  
  // Hyperparameters for the batch scale
  rho = _rho;
  theta = _theta;

  // Set the size of the objects to hold the component specific parameters
  mu.set_size(P, K);
  mu.zeros();
  
  cov.set_size(P, P, K);
  cov.zeros();
  
  // Set the size of the objects to hold the batch specific parameters
  m.set_size(P, B);
  m.zeros();
  
  // We are assuming a diagonal structure in the batch scale
  S.set_size(P, B);
  S.zeros();
  
  // Count the number of times proposed values are accepted
  cov_count = zeros<uvec>(K);
  mu_count = zeros<uvec>(K);
  m_count = zeros<uvec>(B);
  S_count = zeros<uvec>(B);
  
  // These will hold vertain matrix operations to avoid computational burden
  // The log determinant of each cluster covariance
  cov_log_det = zeros<vec>(K);
  
  // The log determinant of the covariance combination
  cov_comb_log_det.set_size(K, B);
  cov_comb_log_det.zeros();
  
  // The possible combinations for the sum of the cluster and batch means
  mean_sum.set_size(P, K * B);
  mean_sum.zeros();
  
  // The combination of each possible cluster and batch covariance
  cov_comb.set_size(P, P, K * B);
  cov_comb.zeros();
  
  // Inverse of the cluster covariance
  cov_inv.set_size(P, P, K);
  cov_inv.zeros();
  
  // The inverse of the covariance combination
  cov_comb_inv.set_size(P, P, K * B);
  cov_comb_inv.zeros();
  
  // The batch corrected data
  Y.set_size(N, P);
  Y.zeros();
  
  // The proposal windows for the cluster and batch parameters
  mu_proposal_window = _mu_proposal_window;
  cov_proposal_window = _cov_proposal_window;
  m_proposal_window = _m_proposal_window;
  S_proposal_window = _S_proposal_window;
  
  // flag indicating if the hyperparameter on the batch shift is sampled or not
  sample_m_scale = _sample_m_scale;
};

void mvnSampler::sampleCovPrior() {
  for(uword k = 0; k < K; k++){
    cov.slice(k) = iwishrnd(scale, nu);
  }
};

void mvnSampler::sampleMuPrior() {
  for(uword k = 0; k < K; k++){
    mu.col(k) = mvnrnd(xi, (1.0/kappa) * cov.slice(k), 1);
  }
};

void mvnSampler::sampleSPrior() {
  for(uword b = 0; b < B; b++){
    for(uword p = 0; p < P; p++){
      S(p, b) = S_loc + 1.0 / randg<double>( distr_param(rho, 1.0 / theta ) );
    }
  }
};

void mvnSampler::sampleMPrior() {
  for(uword b = 0; b < B; b++){
    for(uword p = 0; p < P; p++){
      m(p, b) = randn<double>() * batch_shift_prior_precision + batch_shift_prior_mean;
    }
  }
};

void mvnSampler::sampleFromPriors() {
  if(sample_m_scale) {
    sampleMScalePrior();
  }
  sampleCovPrior();
  sampleMuPrior();
  sampleSPrior();
  sampleMPrior();
  if(include_interaction) {
    sampleGammaPrior();
  }
};

void mvnSampler::sampleMScalePrior() {
  lambda_2 = rInvGamma(a, b);
  batch_shift_prior_precision = 1.0 / (delta_2 * lambda_2);
}

void mvnSampler::sampleMScalePosterior() {
  double a_pos = 0.0, b_pos = 0.0;
  a_pos = a + 0.5 * P * B;
  b_pos = 0.5 * accu(pow(m, 2.0)) / (2.0 * delta_2) + b;
  lambda_2 = rInvGamma(a_pos, b_pos);
  batch_shift_prior_precision = 1.0 / (delta_2 * lambda_2);
}

// Update the common matrix manipulations to avoid recalculating N times
void mvnSampler::matrixCombinations() {
  
  for(uword k = 0; k < K; k++) {
    cov_inv.slice(k) = inv_sympd(cov.slice(k));
    cov_log_det(k) = log_det(cov.slice(k)).real();
    for(uword b = 0; b < B; b++) {
      cov_comb.slice(k * B + b) = cov.slice(k);
      for(uword p = 0; p < P; p++) {
        cov_comb.slice(k * B + b)(p, p) *= S(p, b);
      }
      cov_comb_log_det(k, b) = log_det(cov_comb.slice(k * B + b)).real();
      cov_comb_inv.slice(k * B + b) = inv_sympd(cov_comb.slice(k * B + b));

      mean_sum.col(k * B + b) = mu.col(k) + m.col(b);
      if(include_interaction) {
        mean_sum.col(k * B + b) += gamma.slice(b).col(k);
      }
    }
  }
};

// The log likelihood of a item belonging to each cluster given the batch label.
arma::vec mvnSampler::itemLogLikelihood(arma::vec x, arma::uword b) {
  
  double exponent = 0.0;
  vec ll(K), dist_to_mean(P), m_b(B);
  ll.zeros();
  dist_to_mean.zeros();
  m_b = m.col(b);
  
  for(uword k = 0; k < K; k++){
    
    // The exponent part of the MVN pdf
    dist_to_mean = x - mean_sum.col(k * B + b);
    exponent = as_scalar(dist_to_mean.t() * cov_comb_inv.slice(k * B + b) * dist_to_mean);
    
    // Normal log likelihood
    ll(k) = -0.5 *(cov_comb_log_det(k, b) + exponent + (double) P * log(2.0 * M_PI));
  }
  
  return(ll);
};

void mvnSampler::calcBIC(){

  // Each occupied component has a weight, a mean vector and a symmetric
  // covariance matrix; each batch has a shift vector and a scale vector.
  // Empty components contribute no parameters, hence K_occ rather than K.
  BIC = 2 * observed_likelihood - (n_param_cluster * K_occ + n_param_batch * B) * std::log(N);

};

double mvnSampler::groupLikelihood(arma::uvec inds,
                               arma::uvec group_inds,
                               arma::vec cov_det,
                               arma::mat mean_sum,
                               arma::cube cov_inv){
  
  uword c = 0;
  double score = 0.0;
  vec dist_from_mean(P);
  
  for (auto& n : inds) {
    c = group_inds(n);
    dist_from_mean = X_t.col(n) - mean_sum.col(c);
    score += as_scalar(cov_det(c) + (dist_from_mean.t() * cov_inv.slice(c) * dist_from_mean));
  }
  return (-0.5 * score);
}


double mvnSampler::mLogKernel(arma::uword b, arma::vec m_b, arma::mat mean_sum) {
  
  double score = 0.0;
  vec dist_from_mean(P);
  dist_from_mean.zeros();
  
  score = groupLikelihood(batch_ind(b),
    labels,
    cov_comb_log_det.col(b),
    mean_sum,
    cov_comb_inv.slices(KB_inds + b)
  );
  
  for(uword p = 0; p < P; p++) {
    score += -0.5 * (batch_shift_prior_precision * std::pow(m_b(p) - batch_shift_prior_mean, 2.0) );
  }

  return score;
};

double mvnSampler::sLogKernel(arma::uword b, 
                          arma::vec S_b, 
                          arma::vec cov_comb_log_det,
                          arma::cube cov_comb_inv) {
  
  double score = 0.0;
  vec dist_from_mean(P);
  dist_from_mean.zeros();
  
  score = groupLikelihood(batch_ind(b),
    labels,
    cov_comb_log_det,
    mean_sum.cols(KB_inds + b),
    cov_comb_inv
  );
  
  for(uword p = 0; p < P; p++) {
    score +=  -((rho + 1) * std::log(S_b(p) - S_loc) + theta / (S_b(p) - S_loc));
  }
  return score;
};

double mvnSampler::muLogKernel(arma::uword k, arma::vec mu_k, arma::mat mean_sum) {
  
  double score = 0.0;
  uvec cluster_ind = find(labels == k);
  vec dist_from_mean(P);
  
  score = groupLikelihood(cluster_ind,
    batch_vec,
    cov_comb_log_det.row(k).t(),
    mean_sum,
    cov_comb_inv.slices(k * B + B_inds)
  );
  
  score += -0.5 * as_scalar(kappa * ((mu_k - xi).t() *  cov_inv.slice(k) * (mu_k - xi)));
  
  return score;
};

double mvnSampler::covLogKernel(arma::uword k, 
                            arma::mat cov_k, 
                            double cov_log_det,
                            arma::mat cov_inv,
                            arma::vec cov_comb_log_det,
                            arma::cube cov_comb_inv) {
  
  // uword b = 0;
  double score = 0.0;
  vec dist_from_mean(P);
  uvec cluster_ind = find(labels == k);
  
  score = groupLikelihood(cluster_ind,
    batch_vec,
    cov_comb_log_det,
    mean_sum.cols(k * B + B_inds),
    cov_comb_inv
  );
  
  score += -0.5 * ( as_scalar((nu + P + 2) * cov_log_det + kappa * ((mu.col(k) - xi).t() * cov_inv * (mu.col(k) - xi)) + trace(scale * cov_inv)) );
  
  return score;
};

void mvnSampler::batchScaleMetropolis() {
  
  bool next = false;
  double u = 0.0, proposed_model_score = 0.0, acceptance_prob = 0.0, current_model_score = 0.0;
  vec S_proposed(P), proposed_cov_comb_log_det(K);
  cube proposed_cov_comb(P, P, K), proposed_cov_comb_inv(P, P, K);
  
  S_proposed.zeros();
  proposed_cov_comb_log_det.zeros();
  proposed_cov_comb.zeros();
  proposed_cov_comb_inv.zeros();
  
  for(uword b = 0; b < B; b++) {
    
    next = false;
    acceptance_prob = 0.0, proposed_model_score = 0.0, current_model_score = 0.0;
    proposed_cov_comb.zeros();
    
    for(uword p = 0; p < P; p++) {
      if((S(p, b) - S_loc) * S_proposal_window < 0.0){ 
        Rcpp::stop("\n\nCurent batch scale equals S_loc");
      }
      if( (1.0 / S_proposal_window) < 0.0) {
        Rcpp::stop("\n\nCurent batch scale proposal window is some how negative?");
      }
      S_proposed(p) = S_loc + randg( distr_param( (S(p, b) - S_loc) * S_proposal_window, 1.0 / S_proposal_window) );
      
      if(S_proposed(p) <= S_loc) {
        next = true;
      }

      // Metropolis-Hastings correction for the asymmetric Gamma proposal.
      // proposed_model_score must accumulate log pi(proposed) + log q(current | proposed)
      // (the REVERSE proposal density: current evaluated under a proposal
      // centred at proposed), and current_model_score must accumulate
      // log pi(current) + log q(proposed | current) (the FORWARD density,
      // i.e. the density of the draw actually taken to get S_proposed).
      // Swapping these (as a previous version of this code did) silently
      // biases the chain: verified by simulating this exact proposal
      // against a known Gamma(5, 2) target, where the swapped assignment
      // recovers a posterior mean of ~1.41 instead of the true 2.5.
      current_model_score += gammaLogLikelihood(S_proposed(p) - S_loc, (S(p, b) - S_loc) * S_proposal_window, S_proposal_window);
      proposed_model_score += gammaLogLikelihood(S(p, b) - S_loc, (S_proposed(p) - S_loc) * S_proposal_window, S_proposal_window);
    }
    
    if(next) {
      continue;
    }
    
    proposed_cov_comb = cov;
    for(uword k = 0; k < K; k++) {
      for(uword p = 0; p < P; p++) {
        proposed_cov_comb.slice(k)(p, p) *= S_proposed(p);
      }
      proposed_cov_comb_log_det(k) = log_det(proposed_cov_comb.slice(k)).real();
      proposed_cov_comb_inv.slice(k) = inv_sympd(proposed_cov_comb.slice(k));
    }
    
    proposed_model_score += sLogKernel(b, 
      S_proposed, 
      proposed_cov_comb_log_det,
      proposed_cov_comb_inv
    );
    
    current_model_score += sLogKernel(b, 
      S.col(b), 
      cov_comb_log_det.col(b),
      cov_comb_inv.slices(KB_inds + b)
    );
    
    u = randu();
    acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));
    
    if(u < acceptance_prob){
      S.col(b) = S_proposed;
      S_count(b)++;
      
      for(uword k = 0; k < K; k++) {
        cov_comb.slice(k * B + b) = proposed_cov_comb.slice(k);
        cov_comb_log_det(k, b) = proposed_cov_comb_log_det(k);
        cov_comb_inv.slice(k * B + b) = proposed_cov_comb_inv.slice(k);
      }
    }
  }
};

void mvnSampler::batchShiftMetorpolis() {
  
  double u = 0.0, 
    proposed_model_score = 0.0, 
    acceptance_prob = 0.0, 
    current_model_score = 0.0;
  
  vec m_proposed(P);
  mat proposed_mean_sum(P, K);
  m_proposed.zeros();
  
  for(arma::uword b = 0; b < B; b++) {
    for(arma::uword p = 0; p < P; p++){
      // The proposal window is now a diagonal matrix of common entries.
      m_proposed(p) = (arma::randn() * m_proposal_window) + m(p, b);
    }
    
    for(arma::uword k = 0; k < K; k++) {
      proposed_mean_sum.col(k) = mu.col(k) + m_proposed;
      if(include_interaction) {
        proposed_mean_sum.col(k) += gamma.slice(b).col(k);
      }
    }

    // We have a symmetric proposal density so the posterior kernel is the only
    // thing we're interested in
    proposed_model_score = mLogKernel(b, m_proposed, proposed_mean_sum);
    current_model_score = mLogKernel(b, m.col(b), mean_sum.cols(KB_inds + b));
    

    u = arma::randu();
    acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));
    
    // If we accept the proposed value, update our count of acceptances, the 
    // batch shift and the mean sum
    if(u < acceptance_prob){
      m.col(b) = m_proposed;
      m_count(b)++;
      
      for(arma::uword k = 0; k < K; k++) {
        mean_sum.col(k * B + b) = proposed_mean_sum.col(k);
      }
    }
  }
};

void mvnSampler::clusterCovarianceMetropolis() {
  
  double u = 0.0, 
    proposed_model_score = 0.0, 
    acceptance_prob = 0.0, 
    current_model_score = 0.0, 
    proposed_cov_log_det = 0.0;
  
  arma::vec proposed_cov_comb_log_det(B);
  arma::mat cov_proposed(P, P), proposed_cov_inv(P, P);
  arma::cube proposed_cov_comb(P, P, B), proposed_cov_comb_inv(P, P, B);
  
  cov_proposed.zeros();
  proposed_cov_inv.zeros();
  proposed_cov_comb_log_det.zeros();
  proposed_cov_comb.zeros();
  proposed_cov_comb_inv.zeros();
  
  for(arma::uword k = 0; k < K ; k++) {
    
    proposed_cov_comb.zeros();
    
    acceptance_prob = 0.0, proposed_model_score = 0.0, current_model_score = 0.0;
    
    // If no items in the class, sample from the prior distribution.
    if(N_k(k) == 0){
      cov_proposed = arma::iwishrnd(scale, nu);
      proposed_cov_inv = arma::inv_sympd(cov_proposed);
      proposed_cov_log_det = arma::log_det(cov_proposed).real();
      for(arma::uword b = 0; b < B; b++) {
        proposed_cov_comb.slice(b) = cov_proposed;
        for(arma::uword p = 0; p < P; p++) {
          proposed_cov_comb.slice(b)(p, p) *= S(p, b);
        }
        proposed_cov_comb_log_det(b) = arma::log_det(proposed_cov_comb.slice(b)).real();
        proposed_cov_comb_inv.slice(b) = arma::inv_sympd(proposed_cov_comb.slice(b));
      }
    } else {
      
      // Proposal density is a Wishart with an expected value of the current 
      // covariance matrix
      cov_proposed = arma::wishrnd(cov.slice(k) / cov_proposal_window, cov_proposal_window);
      
      // Log probability under the proposal density
      proposed_model_score = wishartLogLikelihood(cov.slice(k), cov_proposed / cov_proposal_window, cov_proposal_window, P);
      current_model_score = wishartLogLikelihood(cov_proposed, cov.slice(k) / cov_proposal_window, cov_proposal_window, P);
      
      proposed_cov_inv = arma::inv_sympd(cov_proposed);
      proposed_cov_log_det = arma::log_det(cov_proposed).real();
      
      for(arma::uword b = 0; b < B; b++) {
        proposed_cov_comb.slice(b) = cov_proposed;
        for(arma::uword p = 0; p < P; p++) {
          proposed_cov_comb.slice(b)(p, p) *= S(p, b);
        }
        proposed_cov_comb_log_det(b) = arma::log_det(proposed_cov_comb.slice(b)).real();
        proposed_cov_comb_inv.slice(b) = arma::inv_sympd(proposed_cov_comb.slice(b));
      }
  
      // The proposed model score is updated by the posterior kernel with the 
      // proposed covariance, it's log determinant, inverse and the new 
      // likelihood covariance
      proposed_model_score += covLogKernel(k, 
        cov_proposed,
        proposed_cov_log_det,
        proposed_cov_inv,
        proposed_cov_comb_log_det,
        proposed_cov_comb_inv
      );
      
      // The current model score is based on the posterior kernel for the current
      // covariance matrix
      current_model_score += covLogKernel(k, 
        cov.slice(k), 
        cov_log_det(k),
        cov_inv.slice(k),
        cov_comb_log_det.row(k).t(),
        cov_comb_inv.slices(k * B + B_inds)
      );
      
      // Accept or reject
      u = arma::randu();
      acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));
      
    }
    if( (u < acceptance_prob) || (N_k(k) == 0) ){
      cov.slice(k) = cov_proposed;
      cov_count(k)++;
      
      cov_inv.slice(k) = proposed_cov_inv;
      cov_log_det(k) = proposed_cov_log_det;
      for(arma::uword b = 0; b < B; b++) {
        cov_comb.slice(k * B + b) = proposed_cov_comb.slice(b);
        cov_comb_log_det(k, b) = proposed_cov_comb_log_det(b);
        cov_comb_inv.slice(k * B + b) = proposed_cov_comb_inv.slice(b);
      }
    }
  }
};

void mvnSampler::clusterMeanMetropolis() {
  
  double u = 0.0, proposed_model_score = 0.0, acceptance_prob = 0.0, current_model_score = 0.0;
  arma::vec mu_proposed(P);
  arma::mat proposed_mean_sum(P, B);
  
  mu_proposed.zeros();
  proposed_mean_sum.zeros();
  
  for(arma::uword k = 0; k < K ; k++) {
    if(N_k(k) == 0){
      mu_proposed = arma::mvnrnd(xi, (1.0/kappa) * cov.slice(k), 1);
      for(arma::uword b = 0; b < B; b++) {
        proposed_mean_sum.col(b) = mu_proposed + m.col(b);
        if(include_interaction) {
          proposed_mean_sum.col(b) += gamma.slice(b).col(k);
        }
      }
    } else {
      for(arma::uword p = 0; p < P; p++){
        // The proposal window is now a diagonal matrix of common entries.
        mu_proposed(p) = (arma::randn() * mu_proposal_window) + mu(p, k);
      }
      for(arma::uword b = 0; b < B; b++) {
        proposed_mean_sum.col(b) = mu_proposed + m.col(b);
        if(include_interaction) {
          proposed_mean_sum.col(b) += gamma.slice(b).col(k);
        }
      }

      // The prior is included in the kernel
      proposed_model_score = muLogKernel(k, mu_proposed, proposed_mean_sum);
      current_model_score = muLogKernel(k, mu.col(k), mean_sum.cols(k * B + B_inds));
      
      u = arma::randu();
      acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));
      
    }
    
    if((u < acceptance_prob) || (N_k(k) == 0)) {
      mu.col(k) = mu_proposed;
      mu_count(k)++;
      
      for(arma::uword b = 0; b < B; b++) {
        mean_sum.col(k * B + b) = proposed_mean_sum.col(b);
      }
      
    }
  }
};

// Block Metropolis-Hastings update for gamma(p, ., .), one feature p at a
// time: propose a whole K x B perturbation, project it onto the
// sum-to-zero subspace (doubleCenterMatrix(), sampler.h/.cpp - keeps the
// proposal, and hence gamma itself, confined to that subspace given it
// starts there), and accept/reject the whole K x B move for that feature
// as one block, scored against the FULL dataset (a single (k,b) cell
// cannot be moved on its own without breaking the sum-to-zero constraint
// for the rest of its row/column, so cell-at-a-time updates are not an
// option here - see the "why" in sampler.h's interaction-term section).
void mvnSampler::interactionMetropolis() {

  arma::uvec combined_group = labels * B + batch_vec;
  arma::uvec all_items = arma::regspace<arma::uvec>(0, N - 1);

  arma::vec cov_det_flat(K * B);
  for(arma::uword k = 0; k < K; k++) {
    for(arma::uword b = 0; b < B; b++) {
      cov_det_flat(k * B + b) = cov_comb_log_det(k, b);
    }
  }

  for(arma::uword p = 0; p < P; p++) {

    arma::mat raw_perturbation(K, B);
    for(arma::uword k = 0; k < K; k++) {
      for(arma::uword b = 0; b < B; b++) {
        raw_perturbation(k, b) = arma::randn() * gamma_proposal_window;
      }
    }
    arma::mat perturbation = doubleCenterMatrix(raw_perturbation);

    arma::mat proposed_mean_sum = mean_sum;
    double prior_current = 0.0, prior_proposed = 0.0;
    arma::mat current_gamma_p(K, B), proposed_gamma_p(K, B);

    for(arma::uword k = 0; k < K; k++) {
      for(arma::uword b = 0; b < B; b++) {
        current_gamma_p(k, b) = gamma(p, k, b);
        proposed_gamma_p(k, b) = current_gamma_p(k, b) + perturbation(k, b);

        prior_current += -0.5 * std::pow(current_gamma_p(k, b), 2.0) / tau2_interaction(p);
        prior_proposed += -0.5 * std::pow(proposed_gamma_p(k, b), 2.0) / tau2_interaction(p);

        proposed_mean_sum(p, k * B + b) = mean_sum(p, k * B + b) - current_gamma_p(k, b) + proposed_gamma_p(k, b);
      }
    }

    // Symmetric (Gaussian, projected) proposal: no Hastings correction needed.
    double current_model_score = groupLikelihood(all_items, combined_group, cov_det_flat, mean_sum, cov_comb_inv) + prior_current;
    double proposed_model_score = groupLikelihood(all_items, combined_group, cov_det_flat, proposed_mean_sum, cov_comb_inv) + prior_proposed;

    double u = arma::randu();
    double acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));

    if(u < acceptance_prob) {
      for(arma::uword k = 0; k < K; k++) {
        for(arma::uword b = 0; b < B; b++) {
          gamma(p, k, b) = proposed_gamma_p(k, b);
        }
      }
      mean_sum.row(p) = proposed_mean_sum.row(p);
      gamma_count(p)++;
    }
  }
};

void mvnSampler::updateBatchCorrectedData() {

  arma::mat mu_mat = mu.cols(labels);
  arma::mat location_correction = m.cols(batch_vec);

  // The interaction term gamma_{k,b} is, by construction (see sampler.h),
  // exactly the part of an item's mean that batch b contributes on top of
  // an additive mu_k + m_b for its specific cluster k - i.e. it is just as
  // much "batch effect" as m_b is, it just happens to depend on cluster
  // identity too. Leaving it in would mean "batch-corrected" data for
  // items in a cluster/batch cell with a real interaction still carries a
  // batch-specific offset. Unlike decomposing gamma against mu/m
  // individually (not identified - see sampler.h), this subtraction is
  // well-defined: gamma_{k,b} for the SPECIFIC (k, b) each item actually
  // belongs to is a single identified quantity once sampled, regardless of
  // how the sum mu_k + m_b + gamma_{k,b} happens to have been split.
  if(include_interaction) {
    for(arma::uword n = 0; n < N; n++) {
      location_correction.col(n) += gamma.slice(batch_vec(n)).col(labels(n));
    }
  }

  Y = ((X_t - mu_mat - location_correction) / sqrt(S.cols(batch_vec)) + mu_mat).t();
}

void mvnSampler::metropolisStep() {

  // Metropolis step for cluster parameters
  clusterCovarianceMetropolis();
  clusterMeanMetropolis();

  // Metropolis step for batch parameters
  if(sample_m_scale) {
    sampleMScalePosterior();
  }
  batchScaleMetropolis();
  batchShiftMetorpolis();

  if(include_interaction) {
    interactionMetropolis();
    sampleTauInteractionPosterior();
  }
};
