// semisupervisedSamplerFactory.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "semisupervisedSamplerFactory.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// semisupervisedSamplerFactory class implementation

std::unique_ptr<semisupervisedSampler> semisupervisedSamplerFactory::createSemisupervisedSampler(samplerType type,
  arma::uword K,
  arma::uword B,
  double mu_proposal_window,
  double cov_proposal_window,
  double m_proposal_window,
  double S_proposal_window,
  double t_df_proposal_window,
  arma::uvec labels,
  arma::uvec batch_vec,
  arma::vec concentration,
  arma::mat X,
  arma::uvec fixed,
  double m_scale,
  double rho,
  double theta,
  bool sample_m_scale
) {
  switch (type) {
  // case G: return std::make_unique<gaussianSampler>(K, labels, concentration, X);

  case MVN: 
    return std::make_unique<mvnPredictive>(K,
      B,
      mu_proposal_window,
      cov_proposal_window,
      m_proposal_window,
      S_proposal_window,
      labels,
      batch_vec,
      concentration,
      X,
      fixed,
      m_scale,
      rho,
      theta,
      sample_m_scale
   );
    
  case MVT: 
    return std::make_unique<mvtPredictive>(K,
      B,
      mu_proposal_window,
      cov_proposal_window,
      m_proposal_window,
      S_proposal_window,
      t_df_proposal_window,
      labels,
      batch_vec,
      concentration,
      X,
      fixed,
      m_scale,
      rho,
      theta,
      sample_m_scale
    );
    
  default: throw "invalid sampler type.";
  }

};