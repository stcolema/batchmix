// predictNewBatchMVNSeparationStrategy.h
// =============================================================================
// include guard
#ifndef PREDICTNEWBATCHMVNSEPARATIONSTRATEGY_H
#define PREDICTNEWBATCHMVNSEPARATIONSTRATEGY_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "mvnSamplerSeparationStrategy.h"

// =============================================================================
// predictNewBatchMVNSeparationStrategy function header

//' @title Predict a new batch (MVN_LKJ), given its own data
//' @description As \code{predictNewBatchMVN()}, but for the LKJ
//' separation-strategy sampler (\code{type = "MVN_LKJ"}) - see that
//' function's documentation for the full design. The cluster covariance
//' is frozen directly via \code{cov_draws} exactly as for MVN; its
//' internal correlation/variance decomposition (R/sigma) is never
//' touched during prediction (\code{rMHStep()}/\code{sigmaMHStep()} are
//' skipped whenever \code{predict_mode} is set - see
//' src/mvnSamplerSeparationStrategy.cpp), since only the resulting
//' covariance matrix, not its decomposition, is needed for prediction.
//' @inheritParams predictNewBatchMVN
// [[Rcpp::export]]
Rcpp::List predictNewBatchMVNSeparationStrategy(
    arma::mat X,
    arma::mat X_new,
    arma::uword K,
    arma::uword B,
    arma::uvec batch_vec,
    arma::uvec fixed_new,
    arma::uvec labels_new_init,
    arma::umat label_draws,
    arma::cube means_draws,
    arma::cube cov_draws,
    arma::cube batch_shift_draws,
    arma::cube batch_scale_draws,
    arma::mat shift_new_init,
    arma::mat scale_new_init,
    double m_scale,
    arma::vec lambda_2_draws,
    bool sample_m_scale,
    double rho,
    double theta,
    arma::vec rho_draws,
    double s_scale_prior_mean,
    bool sample_s_scale,
    arma::uword weight_prior_type,
    arma::mat weights_draws,
    arma::cube eta_logit_init_draws,
    arma::mat gp_beta_draws,
    arma::mat pp_mu_draws,
    arma::mat pp_tau2_draws,
    arma::vec gp_tau2_draws,
    arma::vec gp_length_scale_draws,
    arma::vec batch_coordinates_new,
    double eta_proposal_window,
    double m_proposal_window,
    double S_proposal_window,
    arma::uword n_pred_iter,
    arma::uword pred_burn,
    arma::uword pred_thin
);

#endif /* PREDICTNEWBATCHMVNSEPARATIONSTRATEGY_H */
