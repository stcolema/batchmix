// samplerVaryingWeights.h
// =============================================================================
// include guard
#ifndef SAMPLERVARYINGWEIGHTS_H
#define SAMPLERVARYINGWEIGHTS_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "pdfs.h"
# include "sampler.h"
# include "semisupervisedSampler.h"
# include "genericFunctions.h"

// =============================================================================
// samplerVaryingWeights class

class samplerVaryingWeights : public virtual semisupervisedSampler {
  
private:
  
public:
  
  bool use_log_norm_proposal = true;
  double a = 1.0, b = 1.0, mass_proposal_window = 0.02;
  
  using semisupervisedSampler::semisupervisedSampler;
  
  samplerVaryingWeights(
    arma::uword _K,
    arma::uword _B,
    arma::uvec _labels,
    arma::uvec _batch_vec,
    arma::vec _concentration,
    arma::mat _X,
    arma::uvec _fixed
  ) ;
  
  virtual ~samplerVaryingWeights() { };
  
  // Functions required of all mixture models
  // Generic functions
  void updateWeights() override;
  void updateAllocation() override;
  
  virtual void sampleConcentration();
  virtual void sampleConcentrationK(arma::uword k);
  
};

#endif /* SAMPLERVARYINGWEIGHTS_H */