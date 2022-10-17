// # include <RcppArmadillo.h>
// # include <math.h>
// # include <string>
// # include <iostream>

# include "sampler.h"

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
    arma::mat _X)
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

    // // Allocation probability matrix (only makes sense in predictive models)
    // alloc.set_size(N, K);
    // alloc.zeros();

    // The indices of the members of each batch in the dataset
    batch_ind.set_size(B);
    for(uword b = 0; b < B; b++) {
      batch_ind(b) = find(batch_vec == b);
    }
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

  double u = 0.0;
  uvec uniqueK;
  vec comp_prob(K);

  complete_likelihood = 0.0;
  observed_likelihood = 0.0;
  
  for(uword n = 0; n < N; n++){

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

    labels(n) = sum(u > cumsum(comp_prob));

    // // The allocation really only makes sen
    // alloc.row(n) = comp_prob.t();
    
    // Update the complete likelihood based on the new labelling
    complete_likelihood += ll(labels(n));
  }

  // // The model log likelihood
  // observed_likelihood = accu(likelihood);

  // Number of occupied components (used in BIC calculation)
  uniqueK = unique(labels);
  K_occ = uniqueK.n_elem;
};
