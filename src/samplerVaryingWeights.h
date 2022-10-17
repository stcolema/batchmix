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
  
  bool use_log_norm_proposal = false, model_1_used = false;
  double concentraion_shape_hyper = 1.0, concentraion_rate_hyper = 1.0, mass_proposal_window = 125,
    beta = 0.05, gamma_proposal_window = 125;
  mat normalised_weights;
  
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
  double weightLogPosteriorKernel(
      uword N_kb,
      double weight,
      double mass,
      vec weights
  );
  void sampleWeight(uword k, uword b);
    
  void updateAllocation() override;
  
  virtual void sampleConcentration();
  virtual void sampleConcentrationK(arma::uword k);
  double concentrationLogPosteriorKernel(
      double mass, 
      vec concentration_vec,
      vec weights
  );
  
};

#endif /* SAMPLERVARYINGWEIGHTS_H */