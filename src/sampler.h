// sampler.h
// =============================================================================
// include guard
#ifndef SAMPLER_H
#define SAMPLER_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>

// =============================================================================
// virtual sampler class


//' @name sampler
//' @title sampler
//' @description The virtual sampler class that is the parent to all specific 
//' implementations of the mixture model.
//' @field new Constructor \itemize{
//' \item Parameter: K - the number of components to model.
//' \item Parameter: B - the number of batches present.
//' \item Parameter: labels - N-vector of unsigned integers denoting initial 
//' clustering of the data .
//' \item Parameter: batch_vec - N-vector of unsigned integers denoting 
//' the observed grouping variable.
//' \item Parameter: concentration - K- vector of the prior hyperparameter for 
//' the class weights
//' \item Parameter: X - an N x P matrix of the observed data to model.
//' }
//' @field updateWeights Update the weights of each component based on current
//' clustering.
//' @field updateAllocation Sample a new clustering.
//' @field sampleFromPrior Sample parameter values from the prior distributions.
//' @field metropolisStep Perform Metropolis-Hastings sampling for the class and  
//' batch parameters.
//' @field calcBIC Calculate the BIC of the model.
//' @field logLikelihood Calculate the likelihood of a given data point in each
//' component. \itemize{
//' \item Parameter: x - a data point.
//' \item Parameter: b - the grouping variable category associated with x.
//' }
class sampler {

private:

public:

  
  arma::uword K, B, N, P, K_occ, accepted = 0;
  
  double observed_likelihood = 0.0,
    BIC = 0.0, 
    complete_likelihood = 0.0;
  
  arma::uvec labels, N_k, batch_vec, N_b, KB_inds, B_inds;
  arma::vec concentration, w, ll, likelihood;
  arma::umat members;
  arma::mat X, X_t; //, alloc;
  arma::field<arma::uvec> batch_ind;

  // Parametrised class
  sampler(
    arma::uword _K,
    arma::uword _B,
    arma::uvec _labels,
    arma::uvec _batch_vec,
    arma::vec _concentration,
    arma::mat _X);

  // Destructor (I haven't found a way of defining destructors outside headers)
  virtual ~sampler() { };

  // Functions required of all mixture models
  // Generic functions
  virtual void updateWeights();
  virtual void updateAllocation();
  
  // Mixture specific functions (therefore virtual at this level)
  virtual void metropolisStep() = 0;
  virtual void sampleFromPriors() = 0;
  virtual void calcBIC() = 0;
  virtual arma::vec itemLogLikelihood(arma::vec x, arma::uword b) = 0;

};

#endif /* SAMPLER_H */