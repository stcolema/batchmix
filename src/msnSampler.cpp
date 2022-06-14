// msnSampler.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "pdfs.h"
# include "sampler.h"
# include "semisupervisedSampler.h"
# include "mvnSampler.h"
# include "msnSampler.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// msnSampler class



msnSampler::msnSampler(                           
  arma::uword _K,
  arma::uword _B,
  double _mu_proposal_window,
  double _cov_proposal_window,
  double _m_proposal_window,
  double _S_proposal_window,
  double _phi_proposal_window,
  arma::uvec _labels,
  arma::uvec _batch_vec,
  arma::vec _concentration,
  arma::mat _X,
  double _m_scale,
  double _rho,
  double _theta
) : sampler(_K,
_B,
_labels,
_batch_vec,
_concentration,
_X),
mvnSampler(                           
  _K,
  _B,
  _mu_proposal_window,
  _cov_proposal_window,
  _m_proposal_window,
  _S_proposal_window,
  _labels,
  _batch_vec,
  _concentration,
  _X,
  _m_scale,
  _rho,
  _theta
) {
  
  // Hyperparameter for the prior on the shape of the skew normal
  omega = 0.5;
  
  // The shape of the skew normal
  phi.set_size(P, K);
  phi.zeros();
  
  // Count the number of times proposed values are accepted
  phi_count = arma::zeros<arma::uvec>(K);
  
  // These will hold vertain matrix operations to avoid computational burden
  // The standard deviations of the data
  cov_comb_inv_diag_sqrt.set_size(P, K * B);
  cov_comb_inv_diag_sqrt.zeros();
  
  // The proposal windows for the cluster and batch parameters
  phi_proposal_window = _phi_proposal_window;
};

void msnSampler::sampleDFPrior() {
  for(uword k = 0; k < K; k++){
    for(uword p = 0; p < P; p++){
      phi(p, k) =  randn<double>() * omega;
    }
  }
};

void msnSampler::sampleFromPriors() {
  sampleCovPrior();
  sampleMuPrior();
  sampleDFPrior();
  sampleSPrior();
  sampleMPrior();
};

// Update the common matrix manipulations to avoid recalculating N times
void msnSampler::matrixCombinations() {
  
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
      cov_comb_inv_diag_sqrt.col(k * B + b) = sqrt(cov_comb_inv.slice(k * B + b).diag());
      mean_sum.col(k * B + b) = mu.col(k) + m.col(b);
    }
  }
};

// The log likelihood of a item belonging to each cluster given the batch label.
arma::vec msnSampler::itemLogLikelihood(arma::vec x, arma::uword b) {
  
  double exponent = 0.0;
  vec ll(K), dist_to_mean(P);
  ll.zeros();
  dist_to_mean.zeros();
  
  for(uword k = 0; k < K; k++){
    
    // The exponent part of the MVN pdf
    dist_to_mean = x - mean_sum.col(k * B + b);
    exponent = as_scalar(dist_to_mean.t() * cov_comb_inv.slice(k * B + b) * dist_to_mean);
    
    // Normal log likelihood
    ll(k) = log(2.0) + -0.5 * (cov_comb_log_det(k, b) + exponent + (double) P * log(2.0 * M_PI)); 
    ll(k) += log(normcdf(as_scalar(phi.col(k).t() * cov_comb_inv_diag_sqrt(k * B + b) * dist_to_mean)));
    
  }
  
  return(ll);
};

void msnSampler::calcBIC(){
  
  // Each component has a weight, a mean and shape vector and a symmetric covariance 
  // matrix. Each batch has a mean and standard deviations vector.
  // arma::uword n_param = (2 * P + P * (P + 1) * 0.5) * K_occ + (2 * P) * B;
  // BIC = n_param * std::log(N) - 2 * observed_likelihood;
  // 
  // arma::uword n_param_cluster = 1 + 2 * P + P * (P + 1) * 0.5;
  // arma::uword n_param_batch = 2 * P;
  
  // BIC = 2 * observed_likelihood;
  
  BIC = 2 * observed_likelihood - (n_param_batch + n_param_batch) * std::log(N);
  
  // for(arma::uword k = 0; k < K; k++) {
  //   BIC -= n_param_cluster * std::log(N_k(k) + 1);
  // }
  // for(arma::uword b = 0; b < B; b++) {
  //   BIC -= n_param_batch * std::log(N_b(b)+ 1);
  // }
  
};

double msnSampler::batchLikelihood(
  arma::uvec batch_inds,
  arma::vec cov_det,
  arma::mat mean_sum,
  arma::cube cov_inv,
  arma::mat cov_comb_inv_diag_sqrt
) {
  
  uword k = 0;
  double score = 0.0;
  vec dist_from_mean(P);
  dist_from_mean.zeros(); 
  
  for (auto& n : batch_inds) {
    k = labels(n);
    dist_from_mean = X_t.col(n) - mean_sum.col(k);
    score -= 0.5 * as_scalar(dist_from_mean.t() * cov_comb_inv.slice(k) * dist_from_mean);
    score += log(normcdf(as_scalar(phi.col(k).t() 
        * diagmat(cov_comb_inv_diag_sqrt.col(k))
        * dist_from_mean)
      )
    );
  }
  return score;
};

double msnSampler::mLogKernel(arma::uword b, arma::vec m_b, arma::mat mean_sum) {
  
  double score = 0.0;
  
  score = batchLikelihood(
    batch_ind(b), 
    cov_comb_log_det.col(b),
    mean_sum,
    cov_comb_inv.slices(KB_inds + b),
    cov_comb_inv_diag_sqrt.cols(KB_inds + b)
  );
  
  for(uword p = 0; p < P; p++) {
    score -= 0.5 * t * std::pow(m_b(p) - delta, 2.0) ;
  }
  
  return score;
};

double msnSampler::sLogKernel(arma::uword b, 
  arma::vec S_b, 
  arma::vec cov_comb_log_det,
  arma::mat cov_comb_inv_diag_sqrt,
  arma::cube cov_comb_inv
) {
  
  double score = 0.0;
  
  score = batchLikelihood(
    batch_ind(b), 
    cov_comb_log_det,
    mean_sum.cols(KB_inds + b),
    cov_comb_inv,
    cov_comb_inv_diag_sqrt
  );
  
  for(uword p = 0; p < P; p++) {
    score -=  (rho + 1) * std::log(S_b(p)) + theta / S_b(p);
  }
  
  return score;
};

double msnSampler::muLogKernel(arma::uword k, arma::vec mu_k, arma::mat mean_sum) {
  
  uword b = 0;
  double score = 0.0;
  uvec cluster_ind = find(labels == k);
  vec dist_from_mean(P);
  
  for (auto& n : cluster_ind) {
    b = batch_vec(n);
    dist_from_mean = X_t.col(n) - mean_sum.col(b);
    score -= 0.5 * as_scalar(dist_from_mean.t() * cov_comb_inv.slice(k * B + b) * dist_from_mean);
    score += log(normcdf(as_scalar(phi.col(k).t() * diagmat(cov_comb_inv_diag_sqrt.col(k * B + b)) * dist_from_mean)));
  }
  score -= 0.5 * as_scalar(kappa * ((mu_k - xi).t() *  cov_inv.slice(k) * (mu_k - xi)));
  
  return score;
};

double msnSampler::covLogKernel(arma::uword k, 
                    arma::mat cov_k, 
                    double cov_log_det,
                    arma::mat cov_inv,
                    arma::vec cov_comb_log_det,
                    arma::mat cov_comb_inv_diag_sqrt,
                    arma::cube cov_comb_inv) {
  
  uword b = 0;
  double score = 0.0;
  uvec cluster_ind = find(labels == k);
  vec dist_from_mean(P);
  mat curr_sum(P, P);
  for (auto& n : cluster_ind) {
    b = batch_vec(n);
    dist_from_mean = X_t.col(n) - mean_sum.col(k * B + b);
    score -= 0.5 * as_scalar(cov_comb_log_det(b) + (dist_from_mean.t() * cov_comb_inv.slice(b) * dist_from_mean));
    score += log(normcdf(as_scalar(phi.col(k).t() * diagmat(cov_comb_inv_diag_sqrt.col(b)) * dist_from_mean)));
  }
  score -= 0.5 * as_scalar((nu + P + 2) * cov_log_det + kappa * ((mu.col(k) - xi).t() * cov_inv * (mu.col(k) - xi)) + trace(scale * cov_inv));
  
  return score;
};


double msnSampler::phiLogKernel(arma::uword k, arma::vec phi_k) {
  
  uword b = 0;
  double score = 0.0;
  uvec cluster_ind = find(labels == k);
  vec dist_from_mean(P);
  
  for (auto& n : cluster_ind) {
    b = batch_vec(n);
    dist_from_mean = X_t.col(n) - mean_sum.col(b);
    score += log(normcdf(as_scalar(phi.col(k).t() * diagmat(cov_comb_inv_diag_sqrt.col(k * B + b)) * dist_from_mean)));
    
    
  }
  
  for(uword p = 0; p < P; p++) {
    score += -0.5 * std::pow(phi_k(p), 2.0) / omega;
  }
  
  // std::cout << "\nPhi score: " << score;
  
  return score;
};

void msnSampler::batchScaleMetropolis() {
  
  double u = 0.0, proposed_model_score = 0.0, acceptance_prob = 0.0, current_model_score = 0.0;
  vec S_proposed(P), proposed_cov_comb_log_det(K);
  mat proposed_cov_comb_inv_diag_sqrt(P, K);
  cube proposed_cov_comb(P, P, K), proposed_cov_comb_inv(P, P, K);
  
  S_proposed.zeros();
  proposed_cov_comb_log_det.zeros();
  proposed_cov_comb_inv_diag_sqrt.zeros();
  proposed_cov_comb.zeros();
  proposed_cov_comb_inv.zeros();
  
  for(uword b = 0; b < B ; b++) {
    
    acceptance_prob = 0.0, proposed_model_score = 0.0, current_model_score = 0.0;
    
    for(uword p = 0; p < P; p++) {
      
      // S_proposed(p) = std::exp(randn() * S_proposal_window + log(S(p, b)));
      // std::cout << "\n\nS:" << S(p,b);
      
      S_proposed(p) = randg( distr_param( S(p, b) * S_proposal_window, 1.0 / S_proposal_window) );
      
      // std::cout << "\nS proposed:" << S_proposed(p);
      
      // Asymmetric proposal density
      proposed_model_score += gammaLogLikelihood(S(p, b), S_proposed(p) * S_proposal_window, S_proposal_window);
      current_model_score += gammaLogLikelihood(S_proposed(p), S(p, b) * S_proposal_window, S_proposal_window);
      
      
    }
    
    proposed_cov_comb = cov;
    for(uword k = 0; k < K; k++) {
      // proposed_batch_cov_comb.slice(k) = cov.slice(k); // + diagmat(S.col(b))
      for(uword p = 0; p < P; p++) {
        proposed_cov_comb.slice(k)(p, p) *= S_proposed(p);
      }
      proposed_cov_comb_log_det(k) = log_det(proposed_cov_comb.slice(k)).real();
      proposed_cov_comb_inv.slice(k) = inv_sympd(proposed_cov_comb.slice(k));
      proposed_cov_comb_inv_diag_sqrt.col(k) = sqrt(proposed_cov_comb_inv.slice(k).diag());
    }
    
    proposed_model_score += sLogKernel(b,
                                       S_proposed,
                                       proposed_cov_comb_log_det,
                                       proposed_cov_comb_inv_diag_sqrt,
                                       proposed_cov_comb_inv
    );
    
    current_model_score += sLogKernel(b,
                                      S.col(b),
                                      cov_comb_log_det.col(b),
                                      cov_comb_inv_diag_sqrt.cols(KB_inds + b),
                                      cov_comb_inv.slices(KB_inds + b)
    );
    
    // std::cout << "\nModel scores, current: " << current_model_score << "\nproposed: " << proposed_model_score;
    
    u = randu();
    acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));
    
    // std::cout << "\n\nProposed S:\n"<< S_proposed << "\n\nCurrent S:\n" << S.col(b) << "\n\nProposed score: " << proposed_model_score << "\nCurrent score: " << current_model_score;
    
    if(u < acceptance_prob){
      S.col(b) = S_proposed;
      S_count(b)++;
      
      for(uword k = 0; k < K; k++) {
        cov_comb.slice(k * B + b) = proposed_cov_comb.slice(k);
        cov_comb_log_det(k, b) = proposed_cov_comb_log_det(k);
        cov_comb_inv.slice(k * B + b) = proposed_cov_comb_inv.slice(k);
        cov_comb_inv_diag_sqrt.col(k * B + b) = proposed_cov_comb_inv_diag_sqrt.col(k);
      }
    }
  }
};


void msnSampler::clusterCovarianceMetropolis() {
  
  double u = 0.0, proposed_model_score = 0.0, acceptance_prob = 0.0, current_model_score = 0.0, proposed_cov_log_det = 0.0;
  vec proposed_cov_comb_log_det(B);
  mat cov_proposed(P, P), proposed_cov_inv(P, P), proposed_cov_comb_inv_diag_sqrt(P, B);
  cube proposed_cov_comb(P, P, B), proposed_cov_comb_inv(P, P, B);
  
  cov_proposed.zeros();
  proposed_cov_inv.zeros();
  proposed_cov_comb_log_det.zeros();
  proposed_cov_comb_inv_diag_sqrt.zeros();
  proposed_cov_comb.zeros();
  proposed_cov_comb_inv.zeros();
  
  for(uword k = 0; k < K ; k++) {
    
    proposed_cov_comb.zeros();
    acceptance_prob = 0.0, proposed_model_score = 0.0, current_model_score = 0.0;
    
    if(N_k(k) == 0){
      cov_proposed = iwishrnd(scale, nu);
      proposed_cov_inv = inv_sympd(cov_proposed);
      proposed_cov_log_det = log_det(cov_proposed).real();
      for(uword b = 0; b < B; b++) {
        proposed_cov_comb.slice(b) = cov_proposed; // + diagmat(S.col(b))
        for(uword p = 0; p < P; p++) {
          proposed_cov_comb.slice(b)(p, p) *= S(p, b);
        }
        proposed_cov_comb_log_det(b) = log_det(proposed_cov_comb.slice(b)).real();
        proposed_cov_comb_inv.slice(b) = inv_sympd(proposed_cov_comb.slice(b));
        proposed_cov_comb_inv_diag_sqrt.col(b) = sqrt(proposed_cov_comb.slice(b).diag());
      }
    } else {
      
      cov_proposed = wishrnd(cov.slice(k) / cov_proposal_window, cov_proposal_window);
      
      // Log probability under the proposal density
      proposed_model_score = wishartLogLikelihood(cov.slice(k), cov_proposed / cov_proposal_window, cov_proposal_window, P);
      current_model_score = wishartLogLikelihood(cov_proposed, cov.slice(k) / cov_proposal_window, cov_proposal_window, P);
      
      proposed_cov_inv = inv_sympd(cov_proposed);
      proposed_cov_log_det = log_det(cov_proposed).real();
      for(uword b = 0; b < B; b++) {
        proposed_cov_comb.slice(b) = cov_proposed; // + diagmat(S.col(b))
        for(uword p = 0; p < P; p++) {
          proposed_cov_comb.slice(b)(p, p) *= S(p, b);
        }
        proposed_cov_comb_log_det(b) = log_det(proposed_cov_comb.slice(b)).real();
        proposed_cov_comb_inv.slice(b) = inv_sympd(proposed_cov_comb.slice(b));
        proposed_cov_comb_inv_diag_sqrt.col(b) = sqrt(proposed_cov_comb.slice(b).diag());
      }
      
      // The boolean variables indicate use of the old manipulated matrix or the
      // proposed.
      proposed_model_score += covLogKernel(k,
                                           cov_proposed,
                                           proposed_cov_log_det,
                                           proposed_cov_inv,
                                           proposed_cov_comb_log_det,
                                           proposed_cov_comb_inv_diag_sqrt,
                                           proposed_cov_comb_inv
      );
      
      current_model_score += covLogKernel(k,
                                          cov.slice(k),
                                          cov_log_det(k),
                                          cov_inv.slice(k),
                                          cov_comb_log_det.row(k).t(),
                                          cov_comb_inv_diag_sqrt.cols(k * B + B_inds),
                                          cov_comb_inv.slices(k * B + B_inds)
      );
      
      // std::cout << "\n\nProposed Cov:\n"<< cov_proposed << "\n\nCurrent cov:\n" << cov.slice(k) << "\n\nProposed score: " << proposed_model_score << "\nCurrent score: " << current_model_score;
      
      
      // Accept or reject
      u = randu();
      acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));
    }
    if((u < acceptance_prob) || (N_k(k) == 0)){
      cov.slice(k) = cov_proposed;
      cov_count(k)++;
      
      cov_inv.slice(k) = proposed_cov_inv;
      cov_log_det(k) = proposed_cov_log_det;
      for(uword b = 0; b < B; b++) {
        cov_comb.slice(k * B + b) = proposed_cov_comb.slice(b);
        cov_comb_log_det(k, b) = proposed_cov_comb_log_det(b);
        cov_comb_inv.slice(k * B + b) = proposed_cov_comb_inv.slice(b);
        cov_comb_inv_diag_sqrt.col(k * B + b) = proposed_cov_comb_inv_diag_sqrt.col(b);
      }
    }
  }
};

void msnSampler::clusterShapeMetropolis() {
  
  double u = 0.0, proposed_model_score = 0.0, acceptance_prob = 0.0, current_model_score = 0.0;
  vec phi_proposed(P);
  phi_proposed.zeros();
  
  for(uword k = 0; k < K ; k++) {
    if(N_k(k) == 0){
      for(uword p = 0; p < P; p++) {
        phi_proposed(p) = randn<double>() * omega;
      }
    } else {
      for(uword p = 0; p < P; p++){
        // The proposal window is now a diagonal matrix of common entries.
        phi_proposed(p) = (randn() * phi_proposal_window) + phi(p, k);
      }
      
      
      // The prior is included in the kernel
      proposed_model_score = phiLogKernel(k, phi_proposed);
      current_model_score = phiLogKernel(k, phi.col(k));
      
      u = randu();
      acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));
      
    }
    
    if((u < acceptance_prob) || (N_k(k) == 0)) {
      phi.col(k) = phi_proposed;
      phi_count(k)++;
    }
  }
  
};

void msnSampler::metropolisStep() {
  
  // Metropolis step for cluster parameters
  clusterCovarianceMetropolis();
  clusterMeanMetropolis();
  
  // Update the shape parameter of the skew normal
  clusterShapeMetropolis();
  
  // Metropolis step for batch parameters if more than 1 batch
  // if(B > 1){
  batchScaleMetropolis();
  batchShiftMetorpolis();
  // }
};
