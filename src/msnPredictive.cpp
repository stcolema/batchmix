// msnSampler.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "pdfs.h"
# include "msnPredictive.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// msnPredictive class

msnPredictive::msnPredictive(
  arma::uword _K,
  arma::uword _B,
  double _mu_proposal_window,
  double _cov_proposal_window,
  double _m_proposal_window,
  double _S_proposal_window,
  double _phi_proposal_window,
  // double _rho,
  // double _theta,
  arma::uvec _labels,
  arma::uvec _batch_vec,
  arma::vec _concentration,
  arma::mat _X,
  arma::uvec _fixed,
  double _m_scale,
  double _rho,
  double _theta
  
) : 
  sampler(_K, _B, _labels, _batch_vec, _concentration, _X),
  mvnSampler(                           
    _K,
    _B,
    _mu_proposal_window,
    _cov_proposal_window,
    _m_proposal_window,
    _S_proposal_window,
    _labels,
    _batch_vec,
    _concentration,
    _X,
    _m_scale,
    _rho,
    _theta
  ),
  msnSampler(_K,
    _B,
    _mu_proposal_window,
    _cov_proposal_window,
    _m_proposal_window,
    _S_proposal_window,
    _phi_proposal_window,
    // _rho,
    // _theta,
    _labels,
    _batch_vec,
    _concentration,
    _X,
    _m_scale,
    _rho,
    _theta
  ),
  semisupervisedSampler(_K, _B, _labels, _batch_vec, _concentration, _X, _fixed)
{
};


