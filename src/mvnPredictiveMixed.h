// mvnPredictiveMixed.h
// =============================================================================
// include guard
#ifndef MVNPREDICTIVEMIXED_H
#define MVNPREDICTIVEMIXED_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "sampler.h"
# include "semisupervisedSampler.h"
# include "mvnSamplerSeparationStrategy.h"
# include "mvnSamplerMixed.h"

// =============================================================================
// mvnPredictiveMixed class

//' @name mvnPredictiveMixed
//' @title Predictive (semi-supervised) mixed-data-type MVN mixture
//' @description The semi-supervised counterpart of mvnSamplerMixed
//' (probit/missing/censored columns via a shared latent Gaussian), following
//' the same diamond-inheritance pattern as mvnPredictive: updateAllocation()
//' respects the ``fixed`` labels via semisupervisedSampler, while data
//' augmentation and every other update are exactly mvnSamplerMixed's.
//'
//' One difference from the fully-continuous predictive samplers matters
//' here: updateLatentData() (inherited from mvnSamplerMixed) overwrites the
//' working X/X_t for every item needing augmentation, INCLUDING fixed
//' (label-known) items whose observed columns are binary, missing or
//' censored. That is intentional and correct - a fixed label constrains
//' which cluster an item belongs to, not what its unobserved latent trait
//' value is, so those items still need their latent Gaussian draw updated
//' every sweep from the (fixed) cluster/batch parameters.
class mvnPredictiveMixed : public mvnSamplerMixed, public semisupervisedSampler {

private:

public:

  using mvnSamplerMixed::mvnSamplerMixed;

  mvnPredictiveMixed(
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
    double _eta,
    arma::uvec _column_type,
    arma::umat _censor_code
  );

  virtual ~mvnPredictiveMixed() { };

};

#endif /* MVNPREDICTIVEMIXED_H */
