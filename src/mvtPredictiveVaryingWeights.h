// mvtPredictiveVaryingWeights.h
// =============================================================================
// include guard
#ifndef MVTPREDICTIVEVARYINGWEIGHTS_H
#define MVTPREDICTIVEVARYINGWEIGHTS_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "pdfs.h"
# include "samplerVaryingWeights.h"
# include "mvtSampler.h"

// =============================================================================
// mvtPredictiveVaryingWeights class

class mvtPredictiveVaryingWeights : public mvtSampler, public samplerVaryingWeights {
  
private:
  
public:
  
  using samplerVaryingWeights::samplerVaryingWeights;
  
  mvtPredictiveVaryingWeights(
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
    double _theta
  ) ;
  
  virtual ~mvtPredictiveVaryingWeights() { };
};

#endif /* MVTPREDICTIVEVARYINGWEIGHTS_H */