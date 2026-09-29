// sampleSemisupervisedMVT.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "sampleSemisupervisedMVT.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// sampleSemisupervisedMVT function implementation

Rcpp::List sampleSemisupervisedMVT (
    arma::mat X,
    arma::uword K,
    arma::uword B,
    arma::uvec labels,
    arma::uvec batch_vec,
    arma::uvec fixed,
    double mu_proposal_window,
    double cov_proposal_window,
    double m_proposal_window,
    double S_proposal_window,
    double t_df_proposal_window,
    arma::uword n_iter,
    arma::uword thin,
    arma::vec concentration,
    double m_scale,
    double rho,
    double theta,
    arma::mat initial_mu,
    arma::cube initial_cov,
    arma::vec initial_df,
    arma::mat initial_m,
    arma::mat initial_S,
    bool mu_initialised,
    bool cov_initialised,
    bool df_initialised,
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
    double pp_mu_prior_sd,
    bool sample_s_scale,
    double a_s,
    double b_s,
    double s_scale_proposal_window
) {

  mvtSampler my_sampler(K,
    B,
    mu_proposal_window,
    cov_proposal_window,
    m_proposal_window,
    S_proposal_window,
    t_df_proposal_window,
    labels,
    batch_vec,
    concentration,
    X,
    fixed,
    m_scale,
    rho,
    theta,
    sample_m_scale,
    sample_s_scale
  );
  my_sampler.a_s = a_s;
  my_sampler.b_s = b_s;
  my_sampler.s_scale_proposal_window = s_scale_proposal_window;

  if(batch_coordinates.n_elem != B) {
    batch_coordinates = arma::regspace<arma::vec>(0, B - 1);
  }
  my_sampler.initialiseInteraction(include_interaction, gamma_proposal_window, a_gamma, b_gamma);
  my_sampler.initialiseBatchWeightPrior(weight_prior_type, batch_coordinates, gp_tau2, gp_length_scale, eta_proposal_window, sample_gp_hyperparameters, gp_hyperparameter_proposal_window, pp_tau2_shape, pp_tau2_rate, pp_mu_prior_sd);

  uword P = X.n_cols, N = X.n_rows, n_saved = std::floor(n_iter / thin);

  // The output matrix
  umat class_record(n_saved, X.n_rows);
  class_record.zeros();

  // We save the BIC at each iteration
  vec BIC_record = zeros<vec>(n_saved),
    observed_likelihood = zeros<vec>(n_saved),
    complete_likelihood = zeros<vec>(n_saved),
    lambda_2_saved = m_scale * ones<vec>(n_saved);

  mat weights_saved(n_saved, K), t_df_saved(n_saved, K);
  weights_saved.zeros();
  t_df_saved.zeros();

  cube mean_sum_saved(P, K * B, n_saved),
    mu_saved(P, K, n_saved),
    m_saved(P, B, n_saved),
    cov_saved(P, K * P, n_saved),
    S_saved(P, B, n_saved),
    cov_comb_saved(P, P * K * B, n_saved),
    alloc(N, K, n_saved),
    batch_corrected_data(N, P, n_saved),
    latent_data(N, P, n_saved),
    gamma_saved(P, K * B, n_saved),
    w_batch_saved(B, K, n_saved),
    eta_alr_saved(B, (K > 0) ? K - 1 : 0, n_saved);

  mu_saved.zeros();
  cov_saved.zeros();
  cov_comb_saved.zeros();
  m_saved.zeros();
  S_saved.zeros();
  alloc.zeros();
  batch_corrected_data.zeros();
  latent_data.zeros();
  gamma_saved.zeros();

  vec rho_saved = zeros<vec>(n_saved),
    gp_tau2_saved = zeros<vec>(n_saved),
    gp_length_scale_saved = zeros<vec>(n_saved);
  arma::mat pp_mu_saved((K > 0) ? K - 1 : 0, n_saved, arma::fill::zeros),
    pp_tau2_saved((K > 0) ? K - 1 : 0, n_saved, arma::fill::zeros),
    gp_beta_saved((K > 0) ? K - 1 : 0, n_saved, arma::fill::zeros);

  uword save_int = 0;

  // Sampler from priors
  my_sampler.sampleFromPriors();

  // Pass initial values if any are given
  if(mu_initialised) {
    my_sampler.mu = initial_mu;
  }
  if(cov_initialised) {
    my_sampler.cov = initial_cov;
  }
  if(df_initialised) {
    my_sampler.t_df = initial_df;
  }
  if(m_initialised) {
    my_sampler.m = initial_m;
  }
  if(S_initialised) {
    my_sampler.S = initial_S;
  }

  my_sampler.matrixCombinations();

  arma::uvec prev_mu_count = my_sampler.mu_count;
  arma::uvec prev_cov_count = my_sampler.cov_count;
  arma::uvec prev_m_count = my_sampler.m_count;
  arma::uvec prev_S_count = my_sampler.S_count;
  arma::uvec prev_t_df_count = my_sampler.t_df_count;
  arma::uvec prev_gamma_count = my_sampler.gamma_count;
  arma::uvec prev_eta_count = my_sampler.eta_count;
  arma::uword prev_gp_hyperparameter_count = my_sampler.gp_hyperparameter_count;
  arma::uword prev_s_scale_count = my_sampler.s_scale_count;

  // Iterate over MCMC moves
  for(uword r = 0; r < n_iter; r++){

    Rcpp::checkUserInterrupt();

    // Complete the data given the current parameters (see mvtSampler);
    // this runs for fixed (label-known) items too - see mvtSampler.h.
    my_sampler.updateLatentData();

    my_sampler.updateWeights();

    // Metropolis step for batch parameters
    my_sampler.metropolisStep();

    if(auto_tune && r < n_burn) {
      double n_adapt = (double) (r + 1);
      my_sampler.mu_proposal_window = robbinsMonroUpdate(my_sampler.mu_proposal_window, arma::mean(arma::conv_to<arma::vec>::from(my_sampler.mu_count - prev_mu_count)), 0.234, n_adapt);
      my_sampler.cov_proposal_window = robbinsMonroUpdate(my_sampler.cov_proposal_window, arma::mean(arma::conv_to<arma::vec>::from(my_sampler.cov_count - prev_cov_count)), 0.234, n_adapt);
      my_sampler.m_proposal_window = robbinsMonroUpdate(my_sampler.m_proposal_window, arma::mean(arma::conv_to<arma::vec>::from(my_sampler.m_count - prev_m_count)), 0.234, n_adapt);
      my_sampler.S_proposal_window = robbinsMonroUpdate(my_sampler.S_proposal_window, arma::mean(arma::conv_to<arma::vec>::from(my_sampler.S_count - prev_S_count)), 0.234, n_adapt);
      my_sampler.t_df_proposal_window = robbinsMonroUpdate(my_sampler.t_df_proposal_window, arma::mean(arma::conv_to<arma::vec>::from(my_sampler.t_df_count - prev_t_df_count)), 0.44, n_adapt);
      if(include_interaction) {
        my_sampler.gamma_proposal_window = robbinsMonroUpdate(my_sampler.gamma_proposal_window, arma::mean(arma::conv_to<arma::vec>::from(arma::vectorise(my_sampler.gamma_count - prev_gamma_count))), 0.234, n_adapt);
      }
      if(weight_prior_type > 0 && K > 1) {
        my_sampler.eta_proposal_window = robbinsMonroUpdate(my_sampler.eta_proposal_window, arma::mean(arma::conv_to<arma::vec>::from(my_sampler.eta_count - prev_eta_count)), 0.234, n_adapt);
        if(weight_prior_type == 2 && sample_gp_hyperparameters) {
          double gp_hyper_rate = (double) (my_sampler.gp_hyperparameter_count - prev_gp_hyperparameter_count);
          my_sampler.gp_hyperparameter_proposal_window = robbinsMonroUpdate(my_sampler.gp_hyperparameter_proposal_window, gp_hyper_rate, 0.234, n_adapt);
        }
      }
      if(sample_s_scale) {
        double s_scale_rate = (double) (my_sampler.s_scale_count - prev_s_scale_count);
        my_sampler.s_scale_proposal_window = robbinsMonroUpdate(my_sampler.s_scale_proposal_window, s_scale_rate, 0.234, n_adapt);
      }
    }
    prev_mu_count = my_sampler.mu_count;
    prev_cov_count = my_sampler.cov_count;
    prev_m_count = my_sampler.m_count;
    prev_S_count = my_sampler.S_count;
    prev_t_df_count = my_sampler.t_df_count;
    prev_gamma_count = my_sampler.gamma_count;
    prev_eta_count = my_sampler.eta_count;
    prev_gp_hyperparameter_count = my_sampler.gp_hyperparameter_count;
    prev_s_scale_count = my_sampler.s_scale_count;

    my_sampler.updateAllocation();

    // Record results
    if((r + 1) % thin == 0){

      // Update the BIC for the current model fit
      my_sampler.calcBIC();
      BIC_record( save_int ) = my_sampler.BIC;
      observed_likelihood( save_int ) = my_sampler.observed_likelihood;
      complete_likelihood( save_int ) = my_sampler.complete_likelihood;

      class_record.row( save_int ) = my_sampler.labels.t();
      alloc.slice( save_int ) = my_sampler.alloc;

      weights_saved.row( save_int ) = my_sampler.w.t();
      mu_saved.slice( save_int ) = my_sampler.mu;
      m_saved.slice( save_int ) = my_sampler.m;
      S_saved.slice( save_int ) = my_sampler.S;
      mean_sum_saved.slice( save_int ) = my_sampler.mean_sum;
      t_df_saved.row( save_int ) = my_sampler.t_df.t();

      lambda_2_saved( save_int ) = my_sampler.lambda_2;

      cov_saved.slice ( save_int ) = reshape(mat(my_sampler.cov.memptr(), my_sampler.cov.n_elem, 1, false), P, P * K);
      cov_comb_saved.slice( save_int) = reshape(mat(my_sampler.cov_comb.memptr(), my_sampler.cov_comb.n_elem, 1, false), P, P * K * B);

      gamma_saved.slice( save_int ) = reshape(mat(my_sampler.gamma.memptr(), my_sampler.gamma.n_elem, 1, false), P, K * B);
      w_batch_saved.slice( save_int ) = my_sampler.w_batch;
      eta_alr_saved.slice( save_int ) = my_sampler.eta_alr;
      rho_saved( save_int ) = my_sampler.rho;
      gp_tau2_saved( save_int ) = my_sampler.gp_tau2;
      gp_length_scale_saved( save_int ) = my_sampler.gp_length_scale;
      pp_mu_saved.col( save_int ) = my_sampler.pp_mu;
      pp_tau2_saved.col( save_int ) = my_sampler.pp_tau2;
      gp_beta_saved.col( save_int ) = my_sampler.gp_beta;

      my_sampler.updateBatchCorrectedData();
      batch_corrected_data.slice( save_int ) =  my_sampler.Y;
      latent_data.slice( save_int ) = my_sampler.X;

      save_int++;
    }
  }

  Rcpp::List out;
  out["samples"] = class_record;
  out["means"] = mu_saved;
  out["covariance"] = cov_saved;
  out["batch_shift"] = m_saved;
  out["batch_scale"] = S_saved;
  out["mean_sum"] = mean_sum_saved;
  out["cov_comb"] = cov_comb_saved;
  out["t_df"] = t_df_saved;
  out["weights"] = weights_saved;
  out["cov_acceptance_rate"] = conv_to< vec >::from(my_sampler.cov_count) / n_iter;
  out["mu_acceptance_rate"] = conv_to< vec >::from(my_sampler.mu_count) / n_iter;
  out["S_acceptance_rate"] = conv_to< vec >::from(my_sampler.S_count) / n_iter;
  out["m_acceptance_rate"] = conv_to< vec >::from(my_sampler.m_count) / n_iter;
  out["t_df_acceptance_rate"] = conv_to< vec >::from(my_sampler.t_df_count) / n_iter;
  out["alloc"] = alloc;
  out["observed_likelihood"] = observed_likelihood;
  out["complete_likelihood"] = complete_likelihood;
  out["BIC"] = BIC_record;
  out["batch_corrected_data"] = batch_corrected_data;
  out["latent_data"] = latent_data;
  out["lambda_2"] = lambda_2_saved;
  out["gamma"] = gamma_saved;
  out["tau2_interaction"] = my_sampler.tau2_interaction;
  out["gamma_acceptance_rate"] = arma::conv_to< arma::vec >::from(arma::vectorise(my_sampler.gamma_count)) / n_iter;
  out["w_batch"] = w_batch_saved;
  out["eta_alr"] = eta_alr_saved;
  out["eta_acceptance_rate"] = arma::conv_to< arma::vec >::from(my_sampler.eta_count) / n_iter;
  out["gp_tau2"] = gp_tau2_saved;
  out["gp_length_scale"] = gp_length_scale_saved;
  out["pp_mu"] = pp_mu_saved;
  out["pp_tau2"] = pp_tau2_saved;
  out["gp_beta"] = gp_beta_saved;
  out["weight_prior_type"] = weight_prior_type;
  out["gp_hyperparameter_acceptance_rate"] = (double) my_sampler.gp_hyperparameter_count / n_iter;
  out["final_mu_proposal_window"] = my_sampler.mu_proposal_window;
  out["final_cov_proposal_window"] = my_sampler.cov_proposal_window;
  out["final_m_proposal_window"] = my_sampler.m_proposal_window;
  out["final_S_proposal_window"] = my_sampler.S_proposal_window;
  out["final_t_df_proposal_window"] = my_sampler.t_df_proposal_window;
  out["final_gamma_proposal_window"] = my_sampler.gamma_proposal_window;
  out["final_eta_proposal_window"] = my_sampler.eta_proposal_window;
  out["final_gp_hyperparameter_proposal_window"] = my_sampler.gp_hyperparameter_proposal_window;
  out["rho"] = rho_saved;
  out["s_scale_prior_mean"] = my_sampler.s_scale_prior_mean;
  out["sample_s_scale"] = sample_s_scale;
  out["s_scale_acceptance_rate"] = (double) my_sampler.s_scale_count / n_iter;
  out["final_s_scale_proposal_window"] = my_sampler.s_scale_proposal_window;
  return out;
};
