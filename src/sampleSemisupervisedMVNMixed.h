// sampleSemisupervisedMVNMixed.h
// =============================================================================
// include guard
#ifndef SAMPLESEMISUPERVISEDMVNMIXED_H
#define SAMPLESEMISUPERVISEDMVNMIXED_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "mvnPredictiveMixed.h"

// =============================================================================
// sampleSemisupervisedMVNMixed function header

//' @title Sample semi-supervised mixed continuous/binary/missing/censored
//' MVN mixture model
//' @description The semi-supervised (``fixed`` labels respected) counterpart
//' of sampleMVNMixed. See sampleMVNMixed() for the model.
//' @param X The data matrix (items in rows); see sampleMVNMixed().
//' @param K The number of components to model.
//' @param B The number of batches to model.
//' @param labels Vector item labels to initialise from.
//' @param batch_vec Observed batch labels.
//' @param fixed Binary vector of the items that are fixed in their initial
//' label.
//' @param column_type P-vector: 0 for continuous, 1 for binary (probit).
//' @param censor_code N x P matrix; 0 = not censored, 1 = left-censored,
//' 2 = right-censored.
//' @param mu_proposal_window,r_proposal_window,sigma_proposal_window,m_proposal_window,S_proposal_window
//' Metropolis-Hastings proposal windows.
//' @param R The number of iterations to run.
//' @param thin The thinning factor.
//' @param concentration K-vector, the prior concentration for the class
//' weights.
//' @param m_scale,rho,theta Hyperparameters for the batch shift/scale
//' priors.
//' @param eta LKJ concentration parameter.
//' @param sample_m_scale Should the batch shift hyperparameter be sampled?
//' @param auto_tune,n_burn,include_interaction,gamma_proposal_window,a_gamma,b_gamma,weight_prior_type,batch_coordinates,gp_tau2,gp_length_scale,eta_proposal_window,sample_gp_hyperparameters,gp_hyperparameter_proposal_window,pp_tau2_shape,pp_tau2_rate,pp_mu_prior_sd
//' Auto-tuning, interaction-term and GP-correlated-weight options; see
//' sampleMVN() for the full description of each.
//' @return A named list of MCMC samples and diagnostics.
// [[Rcpp::export]]
Rcpp::List sampleSemisupervisedMVNMixed(
    arma::mat X,
    arma::uword K,
    arma::uword B,
    arma::uvec labels,
    arma::uvec batch_vec,
    arma::uvec fixed,
    arma::uvec column_type,
    arma::umat censor_code,
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
    double eta,
    bool sample_m_scale,
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
);

#endif /* SAMPLESEMISUPERVISEDMVNMIXED_H */
