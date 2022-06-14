// semisupervisedSampler.cpp
// =============================================================================
// included dependencies
# include "sampler.h"
# include "semisupervisedSampler.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// virtual semisupervisedSampler class

semisupervisedSampler::semisupervisedSampler(
  arma::uword _K,
  arma::uword _B,
  arma::uvec _labels,
  arma::uvec _batch_vec,
  arma::vec _concentration,
  arma::mat _X,
  arma::uvec _fixed
) : 
  sampler(_K, _B, _labels, _batch_vec, _concentration, _X) {
  
  arma::uvec fixed_ind(N);
  
  fixed = _fixed;
  N_fixed = arma::sum(fixed);
  fixed_ind = arma::find(_fixed == 1);
  unfixed_ind = find(fixed == 0);
  
  alloc.set_size(N, K);
  alloc.zeros();
  
  for (auto& n : fixed_ind) {
    alloc(n, labels(n)) = 1.0;
  }
};

void semisupervisedSampler::updateAllocation() {
  
  double u = 0.0;
  arma::uvec uniqueK;
  arma::vec comp_prob(K);
  
  // The model likelihoods
  complete_likelihood = 0.0;
  observed_likelihood = 0.0;
  
  // for (auto& n : unfixed_ind) {
  for (uword n = 0; n < N; n++) {
    
    // The mixture-specific log likelihood for each observation in each class
    ll = itemLogLikelihood(X_t.col(n), batch_vec(n));
    
    // Update with weights
    comp_prob = ll + log(w);
    
    // Record the likelihood - this is used to calculate the observed likelihood
    // likelihood(n) = accu(comp_prob);
    observed_likelihood += accu(comp_prob);
    
    // Handle overflow problems and then normalise to convert to probabilities
    comp_prob = exp(comp_prob - max(comp_prob));
    comp_prob = comp_prob / sum(comp_prob);
    
    // Prediction and update
    u = randu<double>( );
    
    if(fixed(n) == 0) {
      labels(n) = sum(u > cumsum(comp_prob));
      
      // The allocation probability for each class
      alloc.row(n) = comp_prob.t();
      
    }
    
    // Update the complete likelihood based on the new labelling
    complete_likelihood += ll(labels(n));
  }
  
  // The observed model log likelihood
  // observed_likelihood = arma::accu(likelihood);
  
  // Number of occupied components (used in BIC calculation)
  uniqueK = arma::unique(labels);
  K_occ = uniqueK.n_elem;
};
