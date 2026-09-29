// mvnSampler.h
// =============================================================================
// include guard
#ifndef MVNSAMPLER_H
#define MVNSAMPLER_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "sampler.h"
# include "genericFunctions.h"

// =============================================================================
// mvnSampler class

// //' @name mvnSampler
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
class mvnSampler: public sampler {
  
public:
  
  bool sample_m_scale = true;

  // Partial pooling for the batch SCALE's own population concentration
  // (opt-in; default false - see sampler.h's predict_mode note and the
  // design comment above sScaleConcentrationMetropolis() in mvnSampler.cpp
  // for why this can only touch the *concentration*, not the *mean*, of
  // the batch-scale prior: a uniform rescaling of every S_b is exactly
  // compensated by an inverse rescaling of every cluster covariance
  // (cov_comb = cov_k * diag(S_b)), the same non-identifiability the
  // batch-shift/cluster-mean pair already has and already fixes by
  // pinning the shift prior's mean at 0 - see batch_shift_prior_mean).
  // s_scale_prior_mean = theta0 / (rho0 - 1) is fixed at construction from
  // the user's ORIGINAL rho/theta (mirroring delta_2's role for the
  // shift); theta is then kept in sync as s_scale_prior_mean * (rho - 1)
  // every time rho moves, so the prior MEAN of S_b - S_loc never changes,
  // only how tightly batches concentrate around it (larger rho = tighter).
  bool sample_s_scale = false;

  arma::uword n_param_cluster = 0,
    n_param_batch = 0;

  double kappa = 0.01,
    nu = 0.0,
    // Hyperparameters for the batch mean
    batch_shift_prior_mean = 0.0,
    batch_shift_prior_precision = 1.0,
    delta_2 = 0.0,
    t = 0.0,
    m_scale = 0.01,
    lambda_2 = 0.01,

    // Hyperparameters for the batch scale. These choices expects sampled
    // values in the range of 1.2 to 2.0 which seems a sensible prior belief.
    // rho becomes a free, estimated parameter when sample_s_scale is true
    // (see above) - theta is then a DEPENDENT quantity, not a second free
    // parameter (see s_scale_prior_mean).
    rho = 3.0,
    theta = 1.0,
    S_loc = 1.0, // this gives the batch scale a support of (1.0, \infty)

    // Hyperparameters for sampling m_scale
    a = 3.0,
    b = 1.0,

    // Hyperparameters for sampling the scale concentration (rho), only
    // used if sample_s_scale is true: (rho - 2) ~ Gamma(a_s, b_s) a priori
    // (rho > 2 throughout, for a finite prior variance of S_b); s_scale_a/b
    // default to a prior mean of rho - 2 = a_s / b_s = 2, i.e. rho = 4,
    // close to this package's own pre-existing rho = 3 default.
    s_scale_prior_mean = 0.0,
    a_s = 2.0,
    b_s = 1.0,
    s_scale_proposal_window = 0.1,

    // Proposal windows (initialised but assigned values by user)
    mu_proposal_window = 0.0,
    cov_proposal_window = 0.0,
    m_proposal_window = 0.0,
    S_proposal_window = 0.0;

  arma::uword s_scale_count = 0;

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
    arma::uvec _fixed,
    double _m_scale,
    double _rho,
    double _theta,
    bool _sample_m_scale,
    bool _sample_s_scale = false
  );
  
  // Destructor
  virtual ~mvnSampler() { };
  
  // Parameter specific priors
  void sampleCovPrior();
  void sampleMuPrior();
  void sampleSPrior();
  void sampleMPrior();
  
  // M_scale hyperparameter
  void sampleMScalePrior();
  void sampleMScalePosterior();

  // Batch-scale concentration hyperparameter (opt-in; see sample_s_scale
  // above). No closed-form Gibbs update exists (unlike lambda_2's conjugate
  // InvGamma-InvGamma pair), so this is a Metropolis-Hastings step.
  virtual void sScaleConcentrationMetropolis();
  
  // Update the common matrix manipulations to avoid recalculating N times
  virtual void matrixCombinations();

  // The Gibbs step that redraws every missing (NaN) entry of X from its
  // full conditional given the other entries of the same item and the
  // current parameters (a systematic-scan Gibbs sweep within each item),
  // writing the result into the inherited X/X_t. Must be called once per
  // MCMC sweep, before updateAllocation()/metropolisStep() (see the driver
  // loop in sampleSemisupervisedMVN.cpp), so every missing entry is
  // redrawn from the current parameter state every sweep - not imputed
  // once and held fixed. See mvnSamplerMixed::updateLatentData() (the
  // Albert & Chib 1993 / Dunson 2000 conditional-Gaussian-via-precision-
  // matrix construction this is a direct port of, binary/censoring
  // branches removed) and mvtSampler::updateLatentData() (the Student-t
  // analogue, overriding this). Items with nothing missing (every item
  // when X has no missing data at all) are skipped via items_to_augment,
  // so this draws nothing - and leaves every existing golden-master test
  // unaffected - whenever X is complete.
  virtual void updateLatentData();
  
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

  // Batch x cluster interaction term (opt-in; see sampler.h). Inert
  // no-ops in effect (gamma stays zero, never called from metropolisStep)
  // unless initialiseInteraction(true, ...) has been called.
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

#endif /* MVNSAMPLER_H */