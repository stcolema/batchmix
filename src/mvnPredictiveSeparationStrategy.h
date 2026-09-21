// mvnPredictiveSeparationStrategy.h
// =============================================================================
// include guard
#ifndef MVNPREDICTIVESEPSTRAT_H
#define MVNPREDICTIVESEPSTRAT_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "sampler.h"
# include "semisupervisedSampler.h"
# include "mvnSamplerSeparationStrategy.h"

// =============================================================================
// mvnPredictiveSeparationStrategy class

//' @name mvnPredictiveSeparationStrategy
//' @title Predictive (semi-supervised) LKJ/separation-strategy MVN mixture
//' @description The semi-supervised counterpart of
//' mvnSamplerSeparationStrategy, following the same diamond-inheritance
//' pattern as mvnPredictive/mvtPredictive: it mixes in semisupervisedSampler
//' so that updateAllocation() respects the ``fixed`` labels, while every
//' other update (LKJ correlation, log-normal marginal scales, batch shift
//' and scale) is exactly mvnSamplerSeparationStrategy's.
class mvnPredictiveSeparationStrategy : public mvnSamplerSeparationStrategy, public semisupervisedSampler {

private:

public:

  using mvnSamplerSeparationStrategy::mvnSamplerSeparationStrategy;

  mvnPredictiveSeparationStrategy(
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
    double _eta
  );

  virtual ~mvnPredictiveSeparationStrategy() { };

};

#endif /* MVNPREDICTIVESEPSTRAT_H */
