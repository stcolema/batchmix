// sampleSemisupervisedMVN.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "sampleSemisupervisedMVN.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// sampleSemisupervisedMVN function implementation

Rcpp::List sampleSemisupervisedMVN (
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
    arma::uword n_iter,
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
) {

  mvnSampler my_sampler(K,
    B,
    mu_proposal_window,
    cov_proposal_window,
    m_proposal_window,
    S_proposal_window,
    labels,
    batch_vec,
    concentration,
    X,
    fixed,
    m_scale,
    rho,
    theta,
    sample_m_scale
  );

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
    lambda_2_saved = zeros<vec>(n_saved);

  mat weights_saved = zeros<mat>(n_saved, K);

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
  gamma_saved.zeros();

  vec gp_tau2_saved = zeros<vec>(n_saved),
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
  arma::uvec prev_gamma_count = my_sampler.gamma_count;
  arma::uvec prev_eta_count = my_sampler.eta_count;
  arma::uword prev_gp_hyperparameter_count = my_sampler.gp_hyperparameter_count;

  // Iterate over MCMC moves
  for(uword r = 0; r < n_iter; r++){

    Rcpp::checkUserInterrupt();

    // Complete the data given the current parameters (see mvnSampler);
    // this runs for fixed (label-known) items too - see mvnSampler.h.
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
    }
    prev_mu_count = my_sampler.mu_count;
    prev_cov_count = my_sampler.cov_count;
    prev_m_count = my_sampler.m_count;
    prev_S_count = my_sampler.S_count;
    prev_gamma_count = my_sampler.gamma_count;
    prev_eta_count = my_sampler.eta_count;
    prev_gp_hyperparameter_count = my_sampler.gp_hyperparameter_count;

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

      lambda_2_saved( save_int ) = my_sampler.lambda_2;

      cov_saved.slice ( save_int ) = reshape(mat(my_sampler.cov.memptr(), my_sampler.cov.n_elem, 1, false), P, P * K);
      cov_comb_saved.slice( save_int) = reshape(mat(my_sampler.cov_comb.memptr(), my_sampler.cov_comb.n_elem, 1, false), P, P * K * B);

      gamma_saved.slice( save_int ) = reshape(mat(my_sampler.gamma.memptr(), my_sampler.gamma.n_elem, 1, false), P, K * B);
      w_batch_saved.slice( save_int ) = my_sampler.w_batch;
      eta_alr_saved.slice( save_int ) = my_sampler.eta_alr;
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

  return(
    List::create(Named("samples") = class_record,
      Named("means") = mu_saved,
      Named("covariance") = cov_saved,
      Named("batch_shift") = m_saved,
      Named("batch_scale") = S_saved,
      Named("mean_sum") = mean_sum_saved,
      Named("cov_comb") = cov_comb_saved,
      Named("weights") = weights_saved,
      Named("cov_acceptance_rate") = conv_to< vec >::from(my_sampler.cov_count) / n_iter,
      Named("mu_acceptance_rate") = conv_to< vec >::from(my_sampler.mu_count) / n_iter,
      Named("S_acceptance_rate") = conv_to< vec >::from(my_sampler.S_count) / n_iter,
      Named("m_acceptance_rate") = conv_to< vec >::from(my_sampler.m_count) / n_iter,
      Named("alloc") = alloc,
      Named("observed_likelihood") = observed_likelihood,
      Named("complete_likelihood") = complete_likelihood,
      Named("BIC") = BIC_record,
      Named("batch_corrected_data") = batch_corrected_data,
      Named("latent_data") = latent_data,
      Named("lambda_2") = lambda_2_saved,
      Named("gamma") = gamma_saved,
      Named("tau2_interaction") = my_sampler.tau2_interaction,
      Named("gamma_acceptance_rate") = arma::conv_to< arma::vec >::from(arma::vectorise(my_sampler.gamma_count)) / n_iter,
      Named("w_batch") = w_batch_saved,
      Named("eta_alr") = eta_alr_saved,
      Named("eta_acceptance_rate") = arma::conv_to< arma::vec >::from(my_sampler.eta_count) / n_iter,
      Named("gp_tau2") = gp_tau2_saved,
      Named("gp_length_scale") = gp_length_scale_saved,
      Named("pp_mu") = pp_mu_saved,
      Named("pp_tau2") = pp_tau2_saved,
      Named("gp_beta") = gp_beta_saved,
      Named("weight_prior_type") = weight_prior_type,
      Named("gp_hyperparameter_acceptance_rate") = (double) my_sampler.gp_hyperparameter_count / n_iter,
      Named("final_mu_proposal_window") = my_sampler.mu_proposal_window,
      Named("final_cov_proposal_window") = my_sampler.cov_proposal_window,
      Named("final_m_proposal_window") = my_sampler.m_proposal_window,
      Named("final_S_proposal_window") = my_sampler.S_proposal_window,
      Named("final_gamma_proposal_window") = my_sampler.gamma_proposal_window,
      Named("final_eta_proposal_window") = my_sampler.eta_proposal_window,
      Named("final_gp_hyperparameter_proposal_window") = my_sampler.gp_hyperparameter_proposal_window
    )
  );

};
