// sampleMVNMixed.h
// =============================================================================
// include guard
#ifndef SAMPLEMVNMIXED_H
#define SAMPLEMVNMIXED_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "mvnSamplerMixed.h"

// =============================================================================
// sampleMVNMixed function header

//' @title Sample a mixture of multivariate normal distributions with batch
//' effects, mixed continuous/binary (probit) columns, and missing/censored
//' data
//' @description Performs MCMC sampling for a mixture model with batch
//' effects, an LKJ prior on the cluster correlation structure, and a
//' shared latent Gaussian layer supporting binary (probit-linked) columns
//' alongside continuous ones, plus missing-at-random and left/right
//' censored continuous entries.
//' @param X The data matrix (items in rows). Continuous columns hold the
//' observed value, or NaN for a missing entry, or the known censoring
//' bound for a censored entry (see censor_code). Binary columns hold 0/1,
//' or NaN for a missing outcome.
//' @param K The number of components to model.
//' @param B The number of batches to model.
//' @param labels N-vector of unsigned integers denoting initial clustering.
//' @param batch_vec N-vector of unsigned integers denoting batch of origin.
//' @param column_type P-vector: 0 for a continuous column, 1 for a binary
//' column observed via a probit link.
//' @param censor_code N x P matrix, only meaningful for continuous
//' columns: 0 = not censored, 1 = left-censored (true value below the
//' recorded X entry), 2 = right-censored (true value above the recorded X
//' entry).
//' @param mu_proposal_window,r_proposal_window,sigma_proposal_window,m_proposal_window,S_proposal_window
//' Metropolis-Hastings proposal windows.
//' @param R The number of iterations to run.
//' @param thin The thinning factor.
//' @param concentration K-vector, the prior concentration for the class
//' weights.
//' @param m_scale,rho,theta Hyperparameters for the batch shift and scale
//' priors.
//' @param eta LKJ concentration parameter for the correlation matrix
//' prior; eta = 1 is uniform over correlation matrices.
//' @param sample_m_scale Should the batch shift hyperparameter be sampled?
//' @return A named list of MCMC samples and diagnostics.
// [[Rcpp::export]]
Rcpp::List sampleMVNMixed(
    arma::mat X,
    arma::uword K,
    arma::uword B,
    arma::uvec labels,
    arma::uvec batch_vec,
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

#endif /* SAMPLEMVNMIXED_H */
