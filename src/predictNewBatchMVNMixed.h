// predictNewBatchMVNMixed.h
// =============================================================================
// include guard
#ifndef PREDICTNEWBATCHMVNMIXED_H
#define PREDICTNEWBATCHMVNMIXED_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "mvnSamplerMixed.h"

// =============================================================================
// predictNewBatchMVNMixed function header

//' @title Predict a new batch (MVN_MIXED), given its own data
//' @description As \code{predictNewBatchMVN()}, but for the mixed-type
//' sampler (\code{type = "MVN_MIXED"}, continuous/binary/censored
//' columns) - see that function's documentation for the full design.
//' \code{column_type} (which columns are continuous/binary/censored) is
//' per-COLUMN metadata fixed at training time and so applies unchanged to
//' \code{X_new}; \code{censor_code_new} is \code{X_new}'s own per-item
//' censoring indicator, in the same convention as the original fit's
//' \code{censor_code}.
//' @inheritParams predictNewBatchMVN
//' @param column_type P-vector as at training time (0 = continuous,
//' 1 = binary/probit, 2 = censored - see \code{sampleSemisupervisedMVNMixed()}).
//' @param censor_code_new N_new x P matrix, \code{X_new}'s own censoring
//' indicator (same convention as the original fit's \code{censor_code}).
// [[Rcpp::export]]
Rcpp::List predictNewBatchMVNMixed(
    arma::mat X,
    arma::mat X_new,
    arma::uword K,
    arma::uword B,
    arma::uvec batch_vec,
    arma::uvec fixed_new,
    arma::uvec labels_new_init,
    arma::uvec column_type,
    arma::umat censor_code,
    arma::umat censor_code_new,
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
    arma::uword gp_kernel_type,
    double pp_mu_prior_sd,
    arma::vec batch_coordinates_new,
    double eta_proposal_window,
    double m_proposal_window,
    double S_proposal_window,
    arma::uword n_pred_iter,
    arma::uword pred_burn,
    arma::uword pred_thin
);

#endif /* PREDICTNEWBATCHMVNMIXED_H */
