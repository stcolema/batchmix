// mvnSamplerSeparationStrategy.h
// =============================================================================
// include guard
#ifndef MVNSAMPLERSEPSTRAT_H
#define MVNSAMPLERSEPSTRAT_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "sampler.h"
# include "genericFunctions.h"

// =============================================================================
// mvnSamplerSeparationStrategy class

// //' @name mvnSamplerSeparationStrategy
// //' @title Multivariate Normal mixture type
// //' @description The sampler for the Multivariate Normal mixture model for batch effects.
// //' @field new Constructor \itemize{
// //' \item Parameter: K - the number of components to model.
// //' \item Parameter: B - the number of batches present.
// //' \item Parameter: labels - N-vector of unsigned integers denoting initial 
// //' clustering of the data .
// //' \item Parameter: batch_vec - N-vector of unsigned integers denoting 
// //' the observed grouping variable.
// //' \item Parameter: concentration - K- vector of the prior hyperparameter for 
// //' the class weights
// //' \item Parameter: X - an N x P matrix of the observed data to model.
// //' }
// //' @field updateWeights Update the weights of each component based on current 
// //' clustering.
// //' @field updateAllocation Sample a new clustering. 
// //' @field sampleFromPrior Sample values for the batch and class parameters from
// //' their prior distributions.
// //' @field calcBIC Calculate the BIC of the model.
// //' @field logLikelihood Calculate the log-likelihood of a given data point in each
// //' component. \itemize{
// //' \item Parameter: x - a data point.
// //' \item Parameter: b - the associated batch label.
// //' }
// //' @field updateBatchCorrectedData Transform the observed dataset based on 
// //' sampled parameter values to a batch-corrected dataset.
// //' @return custom mvnSampler class.
class mvnSamplerSeparationStrategy: public sampler {
  
public:
  
  bool sample_m_scale = true;
  
  // 1 (weight) + P (mean) + P*(P+1)/2 (covariance, via the R/sigma
  // decomposition: P*(P-1)/2 free off-diagonal entries of R plus P
  // entries of sigma = P*(P+1)/2, the same count as a general covariance
  // matrix, since (R, sigma) is a bijective reparameterisation of it).
  arma::uword n_param_cluster = 1 + P + P * (P + 1) * 0.5,
    n_param_batch = 2 * P;
  
  // Prior hyperparameters and proposal parameters
  double kappa = 0.01,

    // LKJ concentration parameter for the correlation matrix prior, R ~
    // LKJ(eta). eta = 1 is uniform over the space of correlation matrices;
    // eta > 1 shrinks correlations towards 0 (Lewandowski, Kurowicka & Joe,
    // 2009).
    eta = 1.0,

    // Hyperparameters for the log-normal prior on the marginal standard
    // deviations, sigma_p ~ LogNormal(beta, sqrt(xi))
    beta = 0.5 * log(0.72),
    xi = 1.0,

    // Hyperparameters for the batch mean
    batch_shift_prior_mean = 0.0,
    batch_shift_prior_precision = 1.0,
    delta_2 = 0.0,
    t = 0.0,
    m_scale = 0.01,
    lambda_2 = 0.01,
    
    // Hyperparameters for the batch scale. These choices expects sampled
    // values in the range of 1.2 to 2.0 which seems a sensible prior belief.
    rho = 3.0,
    theta = 1.0, 
    S_loc = 1.0, // this gives the batch scale a support of (1.0, \infty)
    
    // Hyperparameters for sampling m_scale
    a = 3.0,
    b = 1.0,
    
    // Proposal windows (initialised but assigned values by user)
    mu_proposal_window = 0.0,
    r_proposal_window = 0.0,
    sigma_proposal_window = 0.0,
    m_proposal_window = 0.0,
    S_proposal_window = 0.0;
  
  arma::uvec mu_count, r_count, sigma_count, m_count, S_count, phi_count, rcond_count;
  arma::vec mu_0, r_log_det, cov_log_det, global_mean, z_k;
  arma::mat scale, mu, m, S, phi, cov_comb_log_det, mean_sum, global_cov, sigma, Y;
  arma::cube R, Sigma_mat, cov, cov_inv, cov_comb, cov_comb_inv;
  
  using sampler::sampler;
  
  mvnSamplerSeparationStrategy(                           
    arma::uword _K,
    arma::uword _B,
    double _mu_proposal_window,
    double _r_proposal_window,
    double _sigma_proposal_window,
    double _m_proposal_window,
    double _S_proposal_window,
    arma::uvec _labels,
    arma::uvec _batch_vec,
    arma::vec _concentration,
    arma::mat _X,
    arma::uvec _fixed,
    double _m_scale,
    double _rho,
    double _theta,
    bool _sample_m_scale,
    double _eta = 1.0
  );
  
  // Destructor
  virtual ~mvnSamplerSeparationStrategy() { };
  
  // Parameter specific priors
  // (sampleCovPrior/sampleSPrior are virtual so that mvnSamplerMixed can
  // override them to fix the scale of binary/probit columns at 1)
  virtual void sampleCovPrior();
  void sampleMuPrior();
  virtual void sampleSPrior();
  void sampleMPrior();
  
  // M_scale hyperparameter
  void sampleMScalePrior();
  void sampleMScalePosterior();
  
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
  virtual double rLogKernel(uword k,
    double r_log_det,
    vec cov_comb_log_det,
    mat cov,
    mat cov_inverse,
    cube cov_comb_inv
  );
  virtual double sigmaLogKernel(uword k,
                                vec cov_comb_log_det,
                                mat sigma_mat,
                                mat cov,
                                mat cov_inverse,
                                cube cov_comb_inv
  );
  
  // Metropolis-Hastings sampling for the batch scale and class covariance
  virtual void batchScaleMetropolis();
  virtual void rMHStep();
  virtual void sigmaMHStep();
  
  // Metropolis sampling for the batch shift and class mean
  virtual void batchShiftMetorpolis();
  virtual void clusterMeanMetropolis();

  // Batch x cluster interaction term (opt-in; see sampler.h).
  virtual void interactionMetropolis();
  
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

#endif /* MVNSAMPLERSEPSTRAT_H */