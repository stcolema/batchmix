// samplerFactory.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "samplerFactory.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// samplerFactory class implementation
std::unique_ptr<sampler> samplerFactory::createSampler(samplerType type,
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
  double m_scale,
  double rho,
  double theta,
  bool sample_m_scale
) {
  switch (type) {
  // case G: return std::make_unique<gaussianSampler>(K, labels, concentration, X);
  
  case MVN: 
    return std::make_unique<mvnSampler>(K,
      B,
      mu_proposal_window,
      cov_proposal_window,
      m_proposal_window,
      S_proposal_window,
      labels,
      batch_vec,
      concentration,
      X,
      m_scale,
      rho,
      theta,
      sample_m_scale
    );
    
  case MVT: 
    return std::make_unique<mvtSampler>(K,
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
      m_scale,
      rho,
      theta,
      sample_m_scale
    );
    
  default: throw "invalid sampler type.";
  }
};
