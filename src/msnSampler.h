// msnSampler.h
// =============================================================================
// include guard
#ifndef MSNSAMPLER_H
#define MSNSAMPLER_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "sampler.h"
# include "mvnSampler.h"

// =============================================================================
// msnSampler class

//' @name msnSampler
//' @title Multivariate Skew Normal mixture
//' @description The sampler for the Multivariate Skew Normal mixture
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
//' }
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
class msnSampler: virtual public mvnSampler {
  
public:
  
  arma::uword n_param_cluster = 1 + 2 * P + P * (P + 1) * 0.5,
    n_param_batch = 2 * P;
  
  double omega, phi_proposal_window;
  arma::uvec phi_count;
  arma::mat cov_comb_inv_diag_sqrt, phi;
  
  using mvnSampler::mvnSampler;
  
  msnSampler(                           
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
  );
  
  // Destructor
  virtual ~msnSampler() { };
  
  void sampleDFPrior();
  
  virtual void sampleFromPriors() override ;
  
  // Update the common matrix manipulations to avoid recalculating N times
  virtual void matrixCombinations() override ;
  
  // The log likelihood of a item belonging to each cluster given the batch label.
  virtual arma::vec itemLogLikelihood(arma::vec x, arma::uword b) override;
  
  virtual void calcBIC() override;
  
  virtual double batchLikelihood(
    arma::uvec batch_inds,
    arma::vec cov_det,
    arma::mat mean_sum,
    arma::cube cov_inv,
    arma::mat cov_comb_inv_diag_sqrt
  );
  
  // virtual double clusterLikelihood();
  
  double mLogKernel(arma::uword b, arma::vec m_b, arma::mat mean_sum) override; 
  double sLogKernel(arma::uword b, 
    arma::vec S_b, 
    arma::vec cov_comb_log_det,
    arma::mat cov_comb_inv_diag_sqrt,
    arma::cube cov_comb_inv
  );

  double muLogKernel(arma::uword k, arma::vec mu_k, arma::mat mean_sum) override;
  
  double covLogKernel(arma::uword k, 
    arma::mat cov_k, 
    double cov_log_det,
    arma::mat cov_inv,
    arma::vec cov_comb_log_det,
    arma::mat cov_comb_inv_diag_sqrt,
    arma::cube cov_comb_inv
  );
  
  double phiLogKernel(arma::uword k, arma::vec phi_k) ;
  
  virtual void batchScaleMetropolis() override;
  virtual void clusterCovarianceMetropolis() override;
  virtual void clusterShapeMetropolis() ;
  
  virtual void metropolisStep() override;
  
};


#endif /* MSNSAMPLER_H */