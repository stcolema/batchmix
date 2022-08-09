
# include "samplerVaryingWeights.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// Parametrised class
samplerVaryingWeights::samplerVaryingWeights(
  arma::uword _K,
  arma::uword _B,
  arma::uvec _labels,
  arma::uvec _batch_vec,
  arma::vec _concentration,
  arma::mat _X,
  arma::uvec _fixed
) :
  semisupervisedSampler::semisupervisedSampler(
    _K,
    _B,
    _labels,
    _batch_vec,
    _concentration,
    _X,
    _fixed
  ),
  sampler(_K, _B, _labels, _batch_vec, _concentration, _X)
{
  
  w.reset();
  w.set_size(K, B);
  w.zeros();
};

// Functions required of all mixture models
// Function to update class / mixture weights
void samplerVaryingWeights::updateWeights(){
  uword N_bk = 0;
  // Used locally as the posterior concentration
  double alpha = 0.0;
  for (uword b = 0; b < B; b++) {
    for (uword k = 0; k < K; k++) {
      // Find how many labels have the value
      members.col(k) = labels == k;
      N_k(k) = sum(members.col(k));
      
      N_bk = accu((labels == k) && (batch_vec == b));
      
      // Update weights by sampling from a Gamma distribution
      alpha  = concentration(k) + N_bk;
      w(k, b) = randg( distr_param(alpha, 1.0) );
    }
    // Convert the cluster weights (previously gamma distributed) to Dirichlet
    // distributed by normalising (if K = 2 this is a Beta)
    w.col(b) = w.col(b) / accu(w.col(b));
  }
};

// Sample the class allocations
void samplerVaryingWeights::updateAllocation() {
  
  uword b = 0;
  double u = 0.0;
  uvec uniqueK;
  vec comp_prob(K);
  
  complete_likelihood = 0.0;
  observed_likelihood = 0.0;
  
  for(uword n = 0; n < N; n++){
    
    b = batch_vec(n);
    
    // The mixture-specific log likelihood for each observation in each class
    ll = itemLogLikelihood(X_t.col(n), b);
    
    // Update with weights
    comp_prob = ll + log(w.col(b));
    
    // Record the likelihood - this is used to calculate the observed likelihood
    // likelihood(n) = accu(comp_prob);
    observed_likelihood += accu(comp_prob);
    
    // Handle overflow problems and then normalise to convert to probabilities
    comp_prob = exp(comp_prob - max(comp_prob));
    comp_prob = comp_prob / sum(comp_prob);
    
    if(fixed(n) == 0) {
      
      // Prediction and update
      u = randu<double>( );
      
      labels(n) = sum(u > cumsum(comp_prob));
      
      // The allocation probability for each class
      alloc.row(n) = comp_prob.t();
    }
    
    // Update the complete likelihood based on the new labelling
    complete_likelihood += ll(labels(n));
  }
  
  // // The model log likelihood
  // observed_likelihood = accu(likelihood);
  
  // Number of occupied components (used in BIC calculation)
  uniqueK = unique(labels);
  K_occ = uniqueK.n_elem;
};

void samplerVaryingWeights::sampleConcentration() {
  for(uword k = 0; k < K; k++) {
    sampleConcentrationK(k);
  }
};

void samplerVaryingWeights::sampleConcentrationK(uword k) {
    bool accepted = false;
    double
      current_mass = 0.0, 
      proposed_mass = 0.0, 
        
      curr_score = 0.0,
      new_score = 0.0,
        
      acceptance_ratio = 0.0;
    
    vec current_weights(B), concentration_doppel = concentration;
    
    current_weights = w.row(k).t();
    current_mass = concentration(k);
    
    for(uword b = 0; b < B; b++) {
      curr_score += lgamma(current_mass) 
        + (current_mass - 1) * log(current_weights(b)) 
        + (a - 1) * log(current_mass) 
        - b * current_mass;
    }
    curr_score += accu(lgamma(concentration));
    
    proposed_mass = proposeNewNonNegativeValue(
      current_mass,
      mass_proposal_window, 
      use_log_norm_proposal
    );
    // current_mass + randn() * mass_proposal_window;
    if(proposed_mass <= 0.0) {
      acceptance_ratio = 0.0;
    } else {
      concentration_doppel(k) = proposed_mass;
      for(uword b = 0; b < B; b++) {
        new_score += lgamma(proposed_mass) 
          + (proposed_mass - 1) * log(current_weights(b)) 
          + (a - 1) * log(proposed_mass) - b * proposed_mass;
      }
      new_score += accu(lgamma(concentration_doppel));
      acceptance_ratio = exp(new_score - curr_score);
    }
    accepted = metropolisAcceptanceStep(acceptance_ratio);
    if(accepted) {
      concentration(k) = proposed_mass;
    }
};
