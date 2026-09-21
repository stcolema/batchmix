// sampleSemisupervisedMVNSeparationStrategy.h
// =============================================================================
// include guard
#ifndef SAMPLESEMISUPERVISEDMVNSEPARATIONSTRATEGY_H
#define SAMPLESEMISUPERVISEDMVNSEPARATIONSTRATEGY_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "mvnPredictiveSeparationStrategy.h"

// =============================================================================
// sampleSemisupervisedMVNSeparationStrategy function header

//' @title Sample semi-supervised LKJ/separation-strategy MVN mixture model
//' @description The semi-supervised (``fixed`` labels respected) counterpart
//' of sampleMVNSeparationStrategy: an LKJ prior on the cluster correlation
//' structure and log-normal marginal scales, rather than the original
//' package's Inverse-Wishart covariance prior.
//' @param X The data matrix to perform clustering upon (items to cluster in
//' rows).
//' @param K The number of components to model.
//' @param B The number of batches to model.
//' @param labels Vector item labels to initialise from.
//' @param batch_vec Observed batch labels.
//' @param fixed Binary vector of the items that are fixed in their initial
//' label.
//' @param mu_proposal_window,r_proposal_window,sigma_proposal_window,m_proposal_window,S_proposal_window
//' Metropolis-Hastings proposal windows.
//' @param R The number of iterations to run for.
//' @param thin thinning factor for samples recorded.
//' @param concentration Vector of concentrations for mixture weights.
//' @param m_scale,rho,theta Hyperparameters for the batch shift/scale priors.
//' @param initial_mu,initial_cov,initial_m,initial_S,mu_initialised,cov_initialised,m_initialised,S_initialised
//' Optional user-supplied initial parameter values.
//' @param sample_m_scale Bool indicating if the batch shift hyperparameter is
//' sampled or fixed.
//' @param eta LKJ concentration parameter; eta = 1 is uniform over
//' correlation matrices.
//' @param auto_tune,n_burn,include_interaction,gamma_proposal_window,a_gamma,b_gamma,weight_prior_type,batch_coordinates,gp_tau2,gp_length_scale,eta_proposal_window,sample_gp_hyperparameters,gp_hyperparameter_proposal_window,pp_tau2_shape,pp_tau2_rate,pp_mu_prior_sd
//' Auto-tuning, interaction-term and GP-correlated-weight options; see
//' sampleMVN() for the full description of each.
//' @return Named list of the different quantities drawn by the sampler.
// [[Rcpp::export]]
Rcpp::List sampleSemisupervisedMVNSeparationStrategy (
    arma::mat X,
    arma::uword K,
    arma::uword B,
    arma::uvec labels,
    arma::uvec batch_vec,
    arma::uvec fixed,
    double mu_proposal_window,
    double r_proposal_window,
    double sigma_proposal_window,
    double m_proposal_window,
    double S_proposal_window,
    arma::uword R,
    arma::uword thin,
    arma::vec concentration,
    double m_scale,
    double rho,
    double theta,
    arma::mat initial_mu,
    arma::cube initial_cov,
    arma::mat initial_m,
    arma::mat initial_S,
    bool mu_initialised,
    bool cov_initialised,
    bool m_initialised,
    bool S_initialised,
    bool sample_m_scale,
    double eta,
    bool auto_tune,
    arma::uword n_burn,
    bool include_interaction,
    double gamma_proposal_window,
    double a_gamma,
    double b_gamma,
    arma::uword weight_prior_type,
    arma::vec batch_coordinates,
    double gp_tau2,
    double gp_length_scale,
    double eta_proposal_window,
    bool sample_gp_hyperparameters,
    double gp_hyperparameter_proposal_window,
    double pp_tau2_shape,
    double pp_tau2_rate,
    double pp_mu_prior_sd
) ;

#endif /* SAMPLESEMISUPERVISEDMVNSEPARATIONSTRATEGY_H */
