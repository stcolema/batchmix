// mvnPredictive.h
// =============================================================================
// include guard
#ifndef MVNPREDITCTIVE_H
#define MVNPREDITCTIVE_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "sampler.h"
# include "semisupervisedSampler.h"
# include "mvnSampler.h"

// =============================================================================
// mvnPredictive class

//' @name mvnPredictive
//' @title Predictive multivariate normal mixture
//' @description The sampler for the semi-supervised Multivariate Normal mixture
//' model for batch effects.
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
//' @field printType Print the sampler type called.
//' @field updateWeights Update the weights of each component based on current 
//' clustering.
//' @field updateAllocation Sample a new clustering. 
//' @field sampleFromPrior Sample values for the batch and class parameters from
//' their prior distributions.
//' @field calcBIC Calculate the BIC of the model.
//' @field logLikelihood Calculate the log-likelihood of a given data point in each
//' component. \itemize{
//' \item Parameter: x - a data point.
//' \item Parameter: b - the associated batch label.
//' }
//' @field updateBatchCorrectedData Transform the observed dataset based on 
//' sampled parameter values to a batch-corrected dataset.
class mvnPredictive : public mvnSampler, public semisupervisedSampler {
  
private:
  
public:
  
  using mvnSampler::mvnSampler;
  
  mvnPredictive(
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
    arma::uvec _fixed,
    double _m_scale,
    double _rho,
    double _theta
  ) ;
  
  virtual ~mvnPredictive() { };
  
};

#endif /* MVNPREDITCTIVE_H */