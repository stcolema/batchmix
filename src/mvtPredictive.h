// mvtPredictive.h
// =============================================================================
// include guard
#ifndef MVTPREDICTIVE_H
#define MVTPREDICTIVE_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "pdfs.h"
# include "sampler.h"
# include "semisupervisedSampler.h"
# include "mvnSampler.h"
# include "mvnPredictive.h"
# include "mvtSampler.h"

// =============================================================================
// mvtPredictive class

class mvtPredictive : public mvtSampler, public semisupervisedSampler {
  
private:
  
public:
  
  using mvtSampler::mvtSampler;
  
  mvtPredictive(
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
    
  virtual ~mvtPredictive() { };
};

#endif /* MVTPREDICTIVE_H */