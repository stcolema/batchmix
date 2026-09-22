// mvtSampler.h
// =============================================================================
// include guard
#ifndef MVTSAMPLER_H
#define MVTSAMPLER_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "sampler.h"
# include "mvnSampler.h"

// =============================================================================
// mvtSampler class


// //' @name mvtSampler
// //' @title Multivariate t mixture sampler
// //' @description The sampler for the Multivariate t mixture
// //' model for batch effects.
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
// //' @return mvtsampler class.
class mvtSampler: public mvnSampler {
  
public:
  
  // arma::uword t_df = 4;
  arma::uword n_param_cluster = 0, n_param_batch = 0;
  
  // t degree of freedom hyperparameters (decision from
  // https://statmodeling.stat.columbia.edu/2015/05/17/do-we-have-any-recommendations-for-priors-for-student_ts-degrees-of-freedom-parameter/)
  // This gives a very wide range and a support of [2.0, infty).
  double psi = 2.0, 
    chi = 0.1, 
    t_loc = 2.0,
    
    // Our proposal window
    t_df_proposal_window = 0.0, 
    
    // A value in the pdf defined by the degrees of freedom which we save to // [[Rcpp::export]]
    // avoid recomputing
    pdf_const = 0.0;
  
  arma::uvec t_df_count;
  arma::vec t_df, pdf_coef;
  
  
  using mvnSampler::mvnSampler;
  
  mvtSampler(                           
    arma::uword _K,
    arma::uword _B,
    double _mu_proposal_window,
    double _cov_proposal_window,
    double _m_proposal_window,
    double _S_proposal_window,
    double _t_df_proposal_window,
    arma::uvec _labels,
    arma::uvec _batch_vec,
    arma::vec _concentration,
    arma::mat _X,
    arma::uvec _fixed,
    double _m_scale,
    double _rho,
    double _theta,
    bool _sample_m_scale
  ) ;
  
  // Destructor
  virtual ~mvtSampler() { };

  // Calculate the leading coefficient of the MVT pdf
  double calcPDFCoef(double t_df);
  
  
  virtual void sampleDFPrior();
  
  virtual void sampleFromPriors() override;
  
  // Update the common matrix manipulations to avoid recalculating N times
  virtual void matrixCombinations() override;

  // Missing-data augmentation for the multivariate-t model. A plain
  // conditional-Gaussian draw (as in mvnSampler::updateLatentData()) is
  // NOT the correct full conditional here, because the joint distribution
  // of X is Student-t, not Gaussian. Instead this uses the standard
  // Gaussian-scale-mixture representation of the multivariate t (Liu &
  // Rubin, 1995, "ML estimation of the t distribution using EM and its
  // extensions", Statistica Sinica 5): X_n | u_n ~ N(mean_sum, cov_comb /
  // u_n), u_n ~ Gamma(df_k / 2, df_k / 2), which marginalises to exactly
  // the t_P(mean_sum, cov_comb, df_k) density this sampler already uses
  // everywhere else (itemLogLikelihood(), dfLogKernel(), etc.). Each sweep,
  // for every item needing augmentation: draw u_n from its full conditional
  // given the item's current complete row (Gamma((df_k + P)/2, (df_k +
  // Mahalanobis^2)/2), the standard scale-mixture posterior), then draw the
  // missing coordinates from the Gaussian conditional implied by
  // N(mean_sum, cov_comb / u_n) - algebraically identical to the Gaussian
  // case's conditional-mean formula (the u_n scaling cancels out of the
  // mean, leaving only the conditional variance scaled by 1/u_n). u_n is a
  // transient, per-sweep auxiliary variable local to this step only - it is
  // not retained state, and every other update (clusterDFMetropolis(),
  // clusterCovarianceMetropolis(), etc.) continues to target the true
  // marginal t-density on the now-complete row, exactly as when X has no
  // missing data.
  virtual void updateLatentData() override;
  
  // The log likelihood of a item belonging to each cluster given the batch label.
  virtual arma::vec itemLogLikelihood(arma::vec x, arma::uword b) override;
  
  virtual void calcBIC() override;
  
  double clusterLikelihood(
    double t_df,
    arma::uvec cluster_ind,
    arma::vec cov_det,
    arma::mat mean_sum,
    arma::cube cov_inv
  ) ;
  
  double batchLikelihood(
    arma::uvec batch_inds,
    arma::uvec labels,
    arma::vec cov_det,
    arma::vec t_df,
    arma::mat mean_sum,
    arma::cube cov_inv
  );
  
  virtual double mLogKernel(arma::uword b, arma::vec m_b, arma::mat mean_sum) override;
  
  virtual double sLogKernel(arma::uword b,
    arma::vec S_b,
    arma::vec cov_comb_log_det,
    arma::cube cov_comb_inv
  ) override;
  
  virtual double muLogKernel(arma::uword k, arma::vec mu_k, arma::mat mean_sum) override;
  
  virtual double covLogKernel(arma::uword k, 
    arma::mat cov_k,
    double cov_log_det,
    arma::mat cov_inv,
    arma::vec cov_comb_log_det,
    arma::cube cov_comb_inv
  ) override;
  
  virtual double dfLogKernel(arma::uword k, 
    double t_df,
    double pdf_coef
  );
  
  virtual void clusterDFMetropolis();

  // Batch x cluster interaction term (opt-in; see sampler.h). Uses the
  // Student-t (not Gaussian) location-shift log-likelihood, since t_df is
  // held fixed while gamma is updated.
  virtual double interactionDataLogLikelihood(arma::mat mean_sum_arg);
  virtual void interactionMetropolis();

  virtual void metropolisStep() override;

};

#endif /* MVTSAMPLER_H */