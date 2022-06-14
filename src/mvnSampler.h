// mvnSampler.h
// =============================================================================
// include guard
#ifndef MVNSAMPLER_H
#define MVNSAMPLER_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "sampler.h"

// =============================================================================
// mvnSampler class

//' @name mvnSampler
//' @title Multivariate Normal mixture type
//' @description The sampler for the Multivariate Normal mixture model for batch effects.
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
class mvnSampler: virtual public sampler {
  
public:
  
  arma::uword n_param_cluster = 1 + P + P * (P + 1) * 0.5, 
    n_param_batch = 2 * P;
  
  // Prior hyperparameters and proposal parameters
  double kappa = 0.01, 
    nu = P + 2, 
    
    // Hyperparameters for the batch mean
    delta = 0.0,
    t = 0.0,
    m_scale = 0.01,
    // lambda,
    
    // Hyperparameters for the batch scale. These choices expects sampled
    // values in the range of 1.2 to 2.0 which seems a sensible prior belief.
    rho = 3.0,
    theta = 1.0, 
    S_loc = 1.0, // this gives the batch scale a support of (1.0, \infty)
    
    // Proposal windows (initialised but assigned values by user)
    mu_proposal_window = 0.0,
    cov_proposal_window = 0.0,
    m_proposal_window = 0.0,
    S_proposal_window = 0.0;
  
  arma::uvec mu_count, cov_count, m_count, S_count, phi_count, rcond_count;
  arma::vec xi, cov_log_det, global_mean;
  arma::mat scale, mu, m, S, phi, cov_comb_log_det, mean_sum, global_cov, Y;
  arma::cube cov, cov_inv, cov_comb, cov_comb_inv;
  
  using sampler::sampler;
  
  mvnSampler(                           
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
    double m_scale,
    double rho,
    double theta
  );
  
  // Destructor
  virtual ~mvnSampler() { };
  
  // Parameter specific priors
  void sampleCovPrior();
  void sampleMuPrior();
  void sampleSPrior();
  void sampleMPrior();
  
  // Update the common matrix manipulations to avoid recalculating N times
  virtual void matrixCombinations();
  
  // The likelihood for a specific batch or class
  virtual double groupLikelihood(arma::uvec inds,
                                 arma::uvec group_inds,
                                 arma::vec cov_det,
                                 arma::mat mean_sum,
                                 arma::cube cov_inv);
  
  // The posterior log kernels for each parameters
  virtual double mLogKernel(arma::uword b, arma::vec m_b, arma::mat mean_sum);
  virtual double sLogKernel(arma::uword b, 
    arma::vec S_b, 
    arma::vec cov_comb_log_det,
    arma::cube cov_comb_inv
  );
  
  virtual double muLogKernel(arma::uword k, arma::vec mu_k, arma::mat mean_sum);
  virtual double covLogKernel(arma::uword k, 
    arma::mat cov_k, 
    double cov_log_det,
    arma::mat cov_inv,
    arma::vec cov_comb_log_det,
    arma::cube cov_comb_inv
  );
  
  // Metropolis-Hastings sampling for the batch scale and class covariance
  virtual void batchScaleMetropolis();
  virtual void clusterCovarianceMetropolis();
  
  // Metropolis sampling for the batch shift and class mean
  virtual void batchShiftMetorpolis();
  virtual void clusterMeanMetropolis();
  
  // Update our inferred, batch-corrected dataset based on the current sampled 
  // values
  virtual void updateBatchCorrectedData();
  
  // Used in determining problems - probably unnecessary now.
  // virtual void checkPositiveDefinite(arma::uword r);
  
  // Mixture specific functions
  virtual void metropolisStep() override;
  virtual void sampleFromPriors() override;
  virtual void calcBIC() override;
  virtual arma::vec itemLogLikelihood(arma::vec x, arma::uword b) override;

};

#endif /* MVNSAMPLER_H */