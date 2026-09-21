// samplerFactory.h
// =============================================================================
// include guard
#ifndef SAMPLEMVN_H
#define SAMPLEMVN_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "mvnSampler.h"

// =============================================================================
// sampleMVN function header

//' @title Sample mixture of multivariate normal distributions with batch effects
//' @description Performs MCMC sampling for a mixture model with batch effects.
//' @param X The data matrix to perform clustering upon (items to cluster in rows).
//' @param K The number of components to model (upper limit on the number of 
//' clusters found).
//' @param B The number of batches to model.
//' @param labels Vector item labels to initialise from.
//' @param batch_vec Observed batch labels.
//' @param mu_proposal_window The standard deviation for the Gaussian proposal
//' density of the cluster means.
//' @param cov_proposal_window The degrees of freedom for the Wishart proposal
//' density of the cluster covariances.
//' @param m_proposal_window The standard deviation for the Gaussian proposal
//' density of the batch mean effects.
//' @param S_proposal_window The rate for the Gamma proposal density of the 
//' batch scale.
//' @param R The number of iterations to run for.
//' @param thin thinning factor for samples recorded.
//' @param concentration Vector of concentrations for mixture weights
//' (recommended to be symmetric).
//' @param m_scale The scale hyperparameter for the batch shift prior 
//' distribution.
//' @param rho The shape of the prior distribution for the batch scale.
//' @param theta The scale of the prior distribution for the batch scale.
//' @param initial_mu A P x K matrix of initial values for the class means.
//' @param initial_cov A P x P x K cube of initial values for the class 
//' covariance matrices.
//' @param initial_m A P x B matrix of initial values for the batch shift 
//' effects.
//' @param initial_S A P x B matrix of initial values for the batch scales.
//' @param mu_initialised Bool indicating if the class means are initialised by
//' the user. If ``false`` then initial values are drawn from the prior 
//' distribution.
//' @param cov_initialised Bool indicating if the class covariance matrices are 
//' initialised by the user. If ``false`` then initial values are drawn from the
//' prior distribution.
//' @param m_initialised Bool indicating if the batch shift effects are 
//' initialised by the user. If ``false`` then initial values are drawn from the
//' prior distribution.
//' @param S_initialised Bool indicating if the batch scales are initialised by 
//' the user. If ``false`` then initial values are drawn from the prior 
//' distribution.
//' @param sample_m_scale Bool indicating if the hyperparameter on the batch
//' shift effect is sampled or given as fixed.
//' @param auto_tune Bool; if true, every proposal window is adapted during
//' the first ``n_burn`` iterations via Robbins-Monro diminishing adaptation
//' (see ``robbinsMonroUpdate()``) instead of staying fixed at the value
//' passed in.
//' @param n_burn Number of iterations treated as burn-in for proposal-window
//' adaptation; ignored if ``auto_tune`` is false. Adaptation is frozen after
//' this many iterations so the post-burn-in chain retains the correct
//' stationary distribution.
//' @param include_interaction Bool; if true, add a batch x cluster
//' interaction term to the mean, gamma_{k,b} ~ N(0, tau2_interaction), with
//' tau2_interaction ~ InvGamma(a_gamma, b_gamma) (partial pooling).
//' @param gamma_proposal_window Proposal window (Gaussian RW SD) for the
//' interaction term; ignored if ``include_interaction`` is false.
//' @param a_gamma,b_gamma Shape/rate of the InvGamma hyperprior on the
//' interaction shrinkage variance; ignored if ``include_interaction`` is
//' false.
//' @param weight_prior_type Integer; 0 = "global" (the default) - a single
//' mixture weight vector shared by every batch, exactly the original
//' behaviour. 1 = "partial pooling" - each batch gets its own weight
//' vector, with the K-1 additive-log-ratio (ALR) coordinates drawn
//' exchangeably around a shared, estimated population mean/variance (no
//' assumed order, distance or covariance structure between batches - see
//' ``pp_tau2_shape``/``pp_tau2_rate``/``pp_mu_prior_sd``). 2 = "gp" - as
//' partial pooling, but the ALR coordinates are instead linked by a
//' Gaussian process over ``batch_coordinates``, for batches with a genuine
//' known ordering in time or space.
//' @param batch_coordinates A B-vector of 1-D coordinates for the batches
//' (e.g. time order); if of length 0, defaults to 0, 1, ..., B - 1. Only
//' used if ``weight_prior_type`` is 2.
//' @param gp_tau2,gp_length_scale GP marginal variance and length scale for
//' the batch-weight kernel; only used if ``weight_prior_type`` is 2.
//' @param eta_proposal_window Proposal window for the ALR-coordinate
//' block Metropolis-Hastings update; used if ``weight_prior_type`` is 1 or 2.
//' @param sample_gp_hyperparameters Bool; if true, ``gp_tau2`` and
//' ``gp_length_scale`` are themselves updated by Metropolis-Hastings rather
//' than held fixed; only used if ``weight_prior_type`` is 2.
//' @param gp_hyperparameter_proposal_window Proposal window for the GP
//' hyperparameter update; only used if ``weight_prior_type`` is 2 and
//' ``sample_gp_hyperparameters`` is true.
//' @param pp_tau2_shape,pp_tau2_rate Shape/rate of the InvGamma hyperprior
//' on each ALR coordinate's population variance tau2_j; only used if
//' ``weight_prior_type`` is 1.
//' @param pp_mu_prior_sd Prior standard deviation for each ALR
//' coordinate's population mean mu_j (prior mu_j ~ N(0, pp_mu_prior_sd^2));
//' only used if ``weight_prior_type`` is 1.
//' @return Named list of the different quantities drawn by the sampler.
// [[Rcpp::export]]
Rcpp::List sampleMVN (
    arma::mat X,
    arma::uword K,
    arma::uword B,
    arma::uvec labels,
    arma::uvec batch_vec,
    double mu_proposal_window,
    double cov_proposal_window,
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

#endif /* SAMPLEMVN_H */