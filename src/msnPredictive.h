// msnPredictive.h
// =============================================================================
// include guard
#ifndef MSNPREDICTIVE_H
#define MSNPREDICTIVE_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "semisupervisedSampler.h"
# include "msnSampler.h"

// =============================================================================
// msnPredictive class



class msnPredictive : public msnSampler, public semisupervisedSampler {
  
private:
  
public:
  
  using msnSampler::msnSampler;
  
  msnPredictive(
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
    arma::uvec _fixed,
    double _m_scale,
    double _rho,
    double _theta
  ) ;
  
  virtual ~msnPredictive() { };
  
};

#endif /* MSNPREDICTIVE_H */