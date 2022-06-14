// semisupervisedSampler.h
// =============================================================================
// include guard
#ifndef SEMISUPERVISEDSAMPLER_H
#define SEMISUPERVISEDSAMPLER_H

// =============================================================================
// included dependencies
# include "sampler.h"
# include <RcppArmadillo.h>

// =============================================================================
// virtual semisupervisedSampler class

//' @name semisupervisedSampler
//' @title Semi-supeervised sampler
//' @description The generic semi-supervised sampler class. Basically includes 
//' an extra argument (fixed) and saves one additional object.
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
//' \item Parameter: fixed - a binary N-vector indicating if the nth label is 
//' obeserved.
//' }
//' @field printType Print the sampler type called. Used to check pointers are 
//' working.
//' @field updateAllocation Sample a new clustering, keeping the fixed labels.
class semisupervisedSampler : public virtual sampler {
private:
  
public:
  
  arma::uword N_fixed = 0;
  arma::uvec fixed, unfixed_ind;
  arma::mat alloc;
  
  using sampler::sampler;
  
  // Parametrised constructor
  semisupervisedSampler(
    arma::uword _K,
    arma::uword _B,
    arma::uvec _labels,
    arma::uvec _batch_vec,
    arma::vec _concentration,
    arma::mat _X,
    arma::uvec _fixed
  );
  
  // Destructor
  virtual ~semisupervisedSampler() { };
  
  // The allocation function is slightly different in the semisupervised sampler.
  virtual void updateAllocation() override;
};

#endif /* SEMISUPERVISEDSAMPLER_H */