// predictNewBatchMVT.h
// =============================================================================
// include guard
#ifndef PREDICTNEWBATCHMVT_H
#define PREDICTNEWBATCHMVT_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "mvtSampler.h"

// =============================================================================
// predictNewBatchMVT function header

//' @title Predict a new batch (MVT), given its own data
//' @description As \code{predictNewBatchMVN()}, but for the multivariate
//' t sampler - see that function's documentation for the full design
//' (composition sampling restricted to one new batch via
//' \code{sampler::predict_mode}). The only addition is \code{t_df_draws}:
//' each cluster's degrees-of-freedom parameter is frozen at its posterior
//' draw's value throughout prediction, exactly like the cluster mean/
//' covariance (\code{mvtSampler::metropolisStep()} skips
//' \code{clusterDFMetropolis()} whenever \code{predict_mode} is set - see
//' src/mvtSampler.cpp).
//' @inheritParams predictNewBatchMVN
//' @param t_df_draws n_draws x K matrix of cluster degrees-of-freedom.
//' @param t_df_proposal_window Inert (t_df is never proposed during
//' prediction); accepted only for constructor-signature parity.
// [[Rcpp::export]]
Rcpp::List predictNewBatchMVT(
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
    arma::mat t_df_draws,
    arma::cube batch_shift_draws,
    arma::cube batch_scale_draws,
    arma::mat shift_new_init,
    arma::mat scale_new_init,
    double m_scale,
    arma::vec lambda_2_draws,
    bool sample_m_scale,
    double rho,
    double theta,
    arma::uword weight_prior_type,
    arma::mat weights_draws,
    arma::cube eta_alr_init_draws,
    arma::mat gp_beta_draws,
    arma::mat pp_mu_draws,
    arma::mat pp_tau2_draws,
    arma::vec gp_tau2_draws,
    arma::vec gp_length_scale_draws,
    arma::vec batch_coordinates_new,
    double eta_proposal_window,
    double m_proposal_window,
    double S_proposal_window,
    double t_df_proposal_window,
    arma::uword n_pred_iter,
    arma::uword pred_burn,
    arma::uword pred_thin
);

#endif /* PREDICTNEWBATCHMVT_H */
