
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
  uword N_kb = 0;
  w.reset();
  w.set_size(K, B);
  w.zeros();
  normalised_weights.set_size(K, B);
  
  for(uword b = 0; b < B; b++) {
    
    for(uword k = 0; k < K; k++) {
      N_kb = accu((labels == k) && (batch_vec == b));
      N_k(k) = accu(labels == k);
      
      concentration(k) = 1.0 + rGamma(
        concentraion_shape_hyper,
        concentraion_rate_hyper
      );
      // w(k, b) = (double) N_kb; // / (double) N_b(b);
      w(k, b) = 1.0 + rGamma(concentration(k), beta);
    }
    // w.fill(10);
    // concentration.fill(10);
    normalised_weights.col(b) = w.col(b) / accu(w.col(b));
  }
};

// Functions required of all mixture models
// Function to update class / mixture weights
void samplerVaryingWeights::updateWeights(){
  // Rcpp::Rcout << "\n\nUpdating weights.";
  uword N_bk = 0;
  // Used locally as the posterior concentration
  double alpha = 0.0;
  for (uword b = 0; b < B; b++) {
    for (uword k = 0; k < K; k++) {
      
      sampleWeight(k, b);
      // Rcpp::Rcout << "\nWeight successfully sampled.";
      
      // Find how many labels have the value
      members.col(k) = labels == k;
      N_k(k) = sum(members.col(k));
      
      // N_bk = accu((labels == k) && (batch_vec == b));
      // 
      // // Update weights by sampling from a Gamma distribution
      // alpha  = (concentration(k) / (double) B) + N_bk;
      // w(k, b) = randg( distr_param(alpha, 1.0) );
    }
    // Convert the cluster weights (previously gamma distributed) to Dirichlet
    // distributed by normalising (if K = 2 this is a Beta)
    // w.col(b) = w.col(b) / accu(w.col(b));
    
    // Rcpp::Rcout << "\nNormalising weights.";
    normalised_weights.col(b) = w.col(b) / accu(w.col(b));
  }
  // normalised_weights = w;
  // Rcpp::Rcout << "\nWeights updated.";
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
    // comp_prob = ll + log(w.col(b));
    comp_prob = ll + log(normalised_weights.col(b));
    
    
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

double samplerVaryingWeights::concentrationLogPosteriorKernel(
    double alpha, 
    vec concentration_vec,
    vec weights
  ) {
  double score = 0.0;
  
  if(model_1_used) {
    score = (concentraion_shape_hyper - 1.0) * log(alpha)
      - concentraion_rate_hyper * alpha
      - (double) B * lgamma(accu(concentration_vec) / (double) B)
      + (double) B * lgamma(alpha / (double) B)
      + ((alpha / (double) B) - 1.0) * accu(log(weights));
  } else {
    score = (concentraion_shape_hyper - 1.0) * log(alpha)
      + alpha * log(beta)
      - (double) B * lgamma(alpha / (double) B)
      - concentraion_rate_hyper * alpha
      + ((alpha / (double) B) - 1.0) * accu(log(weights));
  }
  return score;
}

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
    
    curr_score = concentrationLogPosteriorKernel(
      current_mass, 
      concentration,
      current_weights
    );
    
    // for(uword b = 0; b < B; b++) {
    //   curr_score += lgamma(current_mass) 
    //     + (current_mass - 1) * log(current_weights(b)) 
    //     + (concentraion_shape_hyper - 1) * log(current_mass) 
    //     - concentraion_rate_hyper * current_mass;
    // }
    // curr_score += accu(lgamma(concentration));
    
    proposed_mass = proposeNewNonNegativeValue(
      current_mass,
      mass_proposal_window, 
      use_log_norm_proposal
    );

    if(proposed_mass <= 0.0) {
      accepted = false;
    } else {
      concentration_doppel(k) = proposed_mass;
      
      new_score = concentrationLogPosteriorKernel(
        proposed_mass, 
        concentration_doppel,
        current_weights
      );
      
      // for(uword b = 0; b < B; b++) {
      //   new_score += lgamma(proposed_mass) 
      //     + (proposed_mass - 1) * log(current_weights(b)) 
      //     + (concentraion_shape_hyper - 1) * log(proposed_mass) - concentraion_rate_hyper * proposed_mass;
      // }
      // new_score += accu(lgamma(concentration_doppel));
      
      acceptance_ratio = exp(new_score - curr_score);
      accepted = metropolisAcceptanceStep(acceptance_ratio);
    }
    
    if(accepted) {
      concentration(k) = proposed_mass;
    }
};



double samplerVaryingWeights::weightLogPosteriorKernel(
    uword N_kb,
    double weight,
    double mass,
    vec weights
) {
  double score = 0.0;
  score += ((double) N_kb + (mass / (double) B) - 1.0) * log(weight) 
    - beta * weight
    - (double) N_kb * log(accu(weights));
    
    return score;
}

void samplerVaryingWeights::sampleWeight(uword k, uword b) {
  // Rcpp::Rcout << "\nSampling weight\nk: " << k << "\nb: " << b;
  bool accepted = false;
  uword N_kb = accu((labels == k) && (batch_vec == b));
  double
    current_weight = w(k, b), 
      proposed_weight = 0.0, 
      
      curr_score = 0.0,
      new_score = 0.0,
      
      acceptance_ratio = 0.0;
  
  vec current_weights = w.row(k).t(), proposed_weights = w.row(k).t();
  if(N_kb > 0) {
    if(model_1_used) {
      // if(((concentration(k) / (double) B) + N_kb) <= 0.0) {
      //   Rcpp::Rcout <<  "\n\nConcentration: " << concentration(k);
      //   Rcpp::Rcout <<  "\nConcentration: " << (concentration(k) / (double) B);
      //   Rcpp::Rcout << "\nN_kb: " << N_kb;
      // }
      proposed_weight = rGamma((concentration(k) / (double) B) + (double) N_kb, 1.0);
      accepted = true;
      } else {
      // Rcpp::Rcout << "\nCurrent score.";
      curr_score = weightLogPosteriorKernel(
        N_kb,
        current_weight,
        concentration(k),
        current_weights
      );
      
      // Rcpp::Rcout <<  "\ncurrent_weight: " << current_weight;
      // Rcpp::Rcout <<  "\ngamma_proposal_window: " << gamma_proposal_window;
      // Rcpp::Rcout << "\nuse_log_norm_proposal: " << use_log_norm_proposal;
  
      
      // Rcpp::Rcout << "\nPropose new value.";
      proposed_weight = proposeNewNonNegativeValue(
        current_weight,
        gamma_proposal_window, 
        use_log_norm_proposal
      );
      
      // Rcpp::Rcout << "\nNew value proposed.";
      // current_mass + randn() * mass_proposal_window;
      if(proposed_weight <= 0.0) {
        acceptance_ratio = 0.0;
        accepted = false;
      } else {
        proposed_weights(b) = proposed_weight;
        
        // Rcpp::Rcout << "\nProposed value's score.";
        new_score = weightLogPosteriorKernel(
          N_kb,
          proposed_weight,
          concentration(k),
          proposed_weights
        );
        
        acceptance_ratio = exp(new_score - curr_score);
        accepted = metropolisAcceptanceStep(acceptance_ratio);
      }
    }
    if(accepted) {
      w(k, b) = proposed_weight;
    }
  } else {
    w(k, b) = 1.0 + rGamma(concentration(k), beta);
  }
};
