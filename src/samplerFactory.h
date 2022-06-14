// samplerFactory.h
// =============================================================================
// include guard
#ifndef SAMPLERFACTORY_H
#define SAMPLERFACTORY_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include <memory>
# include "mvnSampler.h"
# include "msnSampler.h"
# include "mvtSampler.h"

// =============================================================================
// samplerFactory class

//' @name samplerFactory
//' @title Factory for different sampler subtypes.
//' @description The factory allows the type of mixture implemented to change 
//' based upon the user input.
//' @field new Constructor \itemize{
//' \item Parameter: samplerType - the density type to be modelled
//' \item Parameter: K - the number of components to model
//' \item Parameter: labels - the initial clustering of the data
//' \item Parameter: concentration - the vector for the prior concentration of 
//' the Dirichlet distribution of the component weights
//' \item Parameter: X - the data to model
//' }
class samplerFactory
{
public:
  enum samplerType {
    // G = 0,
    MVN = 1,
    MVT = 2,
    MSN = 3
  };
  
  static std::unique_ptr<sampler> createSampler(samplerType type,
    arma::uword K,
    arma::uword B,
    double mu_proposal_window,
    double cov_proposal_window,
    double m_proposal_window,
    double S_proposal_window,
    double t_df_proposal_window,
    double phi_proposal_window,
    arma::uvec labels,
    arma::uvec batch_vec,
    arma::vec concentration,
    arma::mat X,
    double m_scale,
    double rho,
    double theta
  );
  
};


#endif /* SAMPLERFACTORY_H */