// predictNewBatchMVNSeparationStrategy.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "predictNewBatchMVNSeparationStrategy.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// predictNewBatchMVNSeparationStrategy function implementation

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
) {

  uword N = X.n_rows, P = X.n_cols, N_new = X_new.n_rows;
  uword B_new = B + 1;
  uword predict_batch = B; // 0-indexed: the new batch is the last one
  uword n_draws = means_draws.n_slices;

  arma::mat X_full = arma::join_cols(X, X_new);

  arma::uvec batch_vec_new_part(N_new);
  batch_vec_new_part.fill(predict_batch);
  arma::uvec batch_vec_full = arma::join_cols(batch_vec, batch_vec_new_part);

  // Every original item is frozen (see sampler.h/predict_mode): marking it
  // fixed = 1 makes updateAllocation() leave its label untouched every
  // sweep, with no other code path needed.
  arma::uvec fixed_full = arma::join_cols(arma::uvec(N, fill::ones), fixed_new);

  arma::uvec init_labels = arma::join_cols(label_draws.col(0), labels_new_init);

  arma::vec concentration = arma::ones<arma::vec>(K); // inert: never resampled during prediction

  mvnSamplerSeparationStrategy my_sampler(
    K, B_new,
    1.0, 1.0, 1.0, m_proposal_window, S_proposal_window, // mu/r/sigma proposal windows are inert (never used) in predict mode
    init_labels, batch_vec_full, concentration, X_full, fixed_full,
    m_scale, rho, theta, sample_m_scale
  );

  my_sampler.initialiseBatchWeightPrior(
    weight_prior_type, batch_coordinates_new,
    gp_tau2_draws.n_elem ? gp_tau2_draws(0) : 1.0,
    gp_length_scale_draws.n_elem ? gp_length_scale_draws(0) : 1.0,
    eta_proposal_window, false, 0.1, 2.0, 1.0, 10.0
  );

  my_sampler.sampleFromPriors();
  my_sampler.predict_mode = true;
  my_sampler.predict_batch = predict_batch;

  uword n_pred_saved = 0;
  for (uword sweep = 1; sweep <= n_pred_iter; sweep++) {
    if (sweep > pred_burn && (sweep - pred_burn) % pred_thin == 0) n_pred_saved++;
  }
  uword n_total = n_draws * n_pred_saved;

  arma::umat label_new_out(n_total, N_new, fill::zeros);
  arma::cube alloc_new_out(N_new, K, n_total, fill::zeros);
  arma::mat shift_new_out(P, n_total, fill::zeros);
  arma::mat scale_new_out(P, n_total, fill::zeros);
  arma::mat weight_new_out(K, n_total, fill::zeros);

  uword save_idx = 0;

  for (uword t = 0; t < n_draws; t++) {

    Rcpp::checkUserInterrupt();

    my_sampler.mu = means_draws.slice(t);
    for (uword k = 0; k < K; k++) {
      my_sampler.cov.slice(k) = cov_draws.slice(t).cols(k * P, (k + 1) * P - 1);
    }

    my_sampler.m.cols(0, B - 1) = batch_shift_draws.slice(t);
    my_sampler.m.col(predict_batch) = shift_new_init.col(t);
    my_sampler.S.cols(0, B - 1) = batch_scale_draws.slice(t);
    my_sampler.S.col(predict_batch) = scale_new_init.col(t);

    my_sampler.lambda_2 = sample_m_scale ? lambda_2_draws(t) : m_scale;
    my_sampler.batch_shift_prior_precision = 1.0 / (my_sampler.delta_2 * my_sampler.lambda_2);

    if (sample_s_scale) {
      my_sampler.rho = rho_draws(t);
      my_sampler.theta = s_scale_prior_mean * (my_sampler.rho - 1.0);
    }

    if (weight_prior_type == 0) {
      my_sampler.w = weights_draws.row(t).t();
    } else {
      my_sampler.eta_logit = eta_logit_init_draws.slice(t);
      if (weight_prior_type == 1) {
        my_sampler.pp_mu = pp_mu_draws.col(t);
        my_sampler.pp_tau2 = pp_tau2_draws.col(t);
      } else {
        my_sampler.gp_beta = gp_beta_draws.col(t);
        my_sampler.gp_tau2 = gp_tau2_draws(t);
        my_sampler.gp_length_scale = gp_length_scale_draws(t);
        double jitter_unused = 0.0;
        my_sampler.buildWellConditionedGPChol(
          my_sampler.gp_tau2, my_sampler.gp_length_scale,
          my_sampler.gp_cov, my_sampler.gp_chol, jitter_unused
        );
      }
      my_sampler.updateSimplexFromLogits();
    }

    my_sampler.matrixCombinations();

    my_sampler.labels = arma::join_cols(label_draws.col(t), labels_new_init);
    for (arma::uword n = N; n < N + N_new; n++) {
      if (fixed_new(n - N) == 1) {
        my_sampler.alloc.row(n).zeros();
        my_sampler.alloc(n, my_sampler.labels(n)) = 1.0;
      }
    }

    for (uword sweep = 1; sweep <= n_pred_iter; sweep++) {

      my_sampler.updateLatentData();
      my_sampler.updateWeights();
      my_sampler.metropolisStep();
      my_sampler.updateAllocation();

      if (sweep > pred_burn && (sweep - pred_burn) % pred_thin == 0) {
        label_new_out.row(save_idx) = my_sampler.labels.subvec(N, N + N_new - 1).t();
        alloc_new_out.slice(save_idx) = my_sampler.alloc.rows(N, N + N_new - 1);
        shift_new_out.col(save_idx) = my_sampler.m.col(predict_batch);
        scale_new_out.col(save_idx) = my_sampler.S.col(predict_batch);
        weight_new_out.col(save_idx) = (weight_prior_type == 0) ?
          my_sampler.w : my_sampler.w_batch.row(predict_batch).t();
        save_idx++;
      }
    }
  }

  Rcpp::List out;
  out["label_new"] = label_new_out;
  out["alloc_new"] = alloc_new_out;
  out["shift_new"] = shift_new_out;
  out["scale_new"] = scale_new_out;
  out["weight_new"] = weight_new_out;
  return out;
};
