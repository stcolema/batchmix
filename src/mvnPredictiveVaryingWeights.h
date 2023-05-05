// mvnPredictiveVaryingWeights.h
// =============================================================================
// include guard
#ifndef MVNPREDICTIVEVARYINGWEIGHTS_H
#define MVNPREDICTIVEVARYINGWEIGHTS_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "pdfs.h"
# include "samplerVaryingWeights.h"
# include "mvnSampler.h"

// =============================================================================
// mvnPredictiveVaryingWeights class

class mvnPredictiveVaryingWeights : public mvnSampler, public samplerVaryingWeights {
  
private:
  
public:
  
  using samplerVaryingWeights::samplerVaryingWeights;
  
  mvnPredictiveVaryingWeights(
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
  
  virtual ~mvnPredictiveVaryingWeights() { };
};

#endif /* MVNPREDICTIVEVARYINGWEIGHTS_H */