// mvtSampler.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "pdfs.h"
# include "sampler.h"
# include "mvnSampler.h"
# include "mvtSampler.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// mvtSampler class

mvtSampler::mvtSampler(                           
  arma::uword _K,
  arma::uword _B,
  double _mu_proposal_window,
  double _cov_proposal_window,
  double _m_proposal_window,
  double _S_proposal_window,
  double _t_df_proposal_window,
  // double _rho,
  // double _theta,
  arma::uvec _labels,
  arma::uvec _batch_vec,
  arma::vec _concentration,
  arma::mat _X,
  double _m_scale,
  double _rho,
  double _theta,
  bool _sample_m_scale
) : sampler(_K,
_B,
_labels,
_batch_vec,
_concentration,
_X), 
mvnSampler(                           
  _K,
  _B,
  _mu_proposal_window,
  _cov_proposal_window,
  _m_proposal_window,
  _S_proposal_window,
  _labels,
  _batch_vec,
  _concentration,
  _X,
  _m_scale,
  _rho,
  _theta,
  _sample_m_scale
) {
  
  n_param_cluster = 2 + P + P * (P + 1) * 0.5;
  n_param_batch = 2 * P;
  
  t_df.set_size(K);
  t_df.zeros();
  
  pdf_coef.set_size(K);
  pdf_coef.zeros();
  
  t_df_count.set_size(K);
  t_df_count.zeros();
  
  // The proposal windows for the cluster and batch parameters
  t_df_proposal_window = _t_df_proposal_window;
};

double mvtSampler::calcPDFCoef(double t_df){
  double x = std::lgamma(0.5 * (t_df + P)) - std::lgamma(0.5 * t_df) - 0.5 * P * log(t_df);
  return x;
};

void mvtSampler::sampleDFPrior() {
  for(arma::uword k = 0; k < K; k++){
    // Draw from a shifted gamma distribution (i.e. gamma with location parameter)
    t_df(k) = t_loc + arma::randg<double>( arma::distr_param(psi, 1.0 / chi));
  }
};

void mvtSampler::sampleFromPriors() {
  if(sample_m_scale) {
    sampleMScalePrior();
  }
  sampleCovPrior();
  sampleMuPrior();
  sampleDFPrior();
  sampleSPrior();
  sampleMPrior();
  if(include_interaction) {
    sampleGammaPrior();
  }
};

// Update the common matrix manipulations to avoid recalculating N times
void mvtSampler::matrixCombinations() {
  
  for(arma::uword k = 0; k < K; k++) {
    pdf_coef(k) = calcPDFCoef(t_df(k));
    cov_inv.slice(k) = arma::inv_sympd(cov.slice(k));
    cov_log_det(k) = arma::log_det(cov.slice(k)).real();
    
    for(arma::uword b = 0; b < B; b++) {
      cov_comb.slice(k * B + b) = cov.slice(k); // + arma::diagmat(S.col(b))
      for(arma::uword p = 0; p < P; p++) {
        cov_comb.slice(k * B + b)(p, p) *= S(p, b);
      }
      cov_comb_log_det(k, b) = arma::log_det(cov_comb.slice(k * B + b)).real();
      cov_comb_inv.slice(k * B + b) = arma::inv_sympd(cov_comb.slice(k * B + b));

      mean_sum.col(k * B + b) = mu.col(k) + m.col(b);
      if(include_interaction) {
        mean_sum.col(k * B + b) += gamma.slice(b).col(k);
      }
    }
  }
};


// The log likelihood of a item belonging to each cluster given the batch label.
arma::vec mvtSampler::itemLogLikelihood(arma::vec x, arma::uword b) {
  
  double u = 0.0, v = 0.0;
  arma::vec ll(K), dist_to_mean(P);
  ll.zeros();
  dist_to_mean.zeros();
  
  for(arma::uword k = 0; k < K; k++){
    
    // The exponent part of the MVN pdf
    dist_to_mean = x - mean_sum.col(k * B + b);
    
    u = arma::as_scalar(dist_to_mean.t() * cov_comb_inv.slice(k * B + b) * dist_to_mean);
    v = (t_df(k) + P) * log(1.0 + (1/t_df(k)) * u);
    
    // MVT Likelihood for cluster k
    ll(k) = pdf_coef(k) - 0.5 * (cov_comb_log_det(k, b) + v + P * log(M_PI));
    
  }
  
  return(ll);
};

void mvtSampler::calcBIC(){

  // Each occupied component has a weight, a mean vector, a symmetric
  // covariance matrix and a degrees-of-freedom parameter; each batch has a
  // shift vector and a scale vector.
  BIC = 2 * observed_likelihood - (n_param_cluster * K_occ + n_param_batch * B) * std::log(N);

};

double mvtSampler::clusterLikelihood(
    double t_df,
    arma::uvec cluster_ind,
    arma::vec cov_det,
    arma::mat mean_sum,
    arma::cube cov_inv
) {
  
  arma::uword b = 0;
  double score = 0.0;
  arma::vec dist_from_mean(P);
  
  for (auto& n : cluster_ind) {
    b = batch_vec(n);
    dist_from_mean = X_t.col(n) - mean_sum.col(b);
    score += cov_det(b) + (t_df + P) * log(1 + (1/t_df) * arma::as_scalar(dist_from_mean.t() * cov_inv.slice(b) * dist_from_mean));
  }
  return (-0.5 * score);
}

double mvtSampler::batchLikelihood(
  arma::uvec batch_inds,
  arma::uvec labels,
  arma::vec cov_det,
  arma::vec t_df,
  arma::mat mean_sum,
  arma::cube cov_inv
) {
  
  arma::uword k = 0;
  double score = 0.0;
  arma::vec dist_from_mean(P);
  
  for (auto& n : batch_inds) {
    k = labels(n);
    dist_from_mean = X_t.col(n) - mean_sum.col(k);
    score += cov_det(k) + (t_df(k) + P) * log(1 + (1/t_df(k)) * arma::as_scalar(dist_from_mean.t() * cov_inv.slice(k) * dist_from_mean));
  }
  return (-0.5 * score);
}

double mvtSampler::mLogKernel(arma::uword b, arma::vec m_b, arma::mat mean_sum) {
  
  double score = 0.0;
  arma::vec dist_from_mean(P);
  dist_from_mean.zeros();
  
  score = batchLikelihood(batch_ind(b), 
    labels, 
    cov_comb_log_det.col(b),
    t_df,
    mean_sum,
    cov_comb_inv.slices(KB_inds + b)
  );
  
  for(arma::uword p = 0; p < P; p++) {
    score += -0.5 * batch_shift_prior_precision * std::pow(m_b(p) - batch_shift_prior_mean, 2.0);
  }
  return score;
};

double mvtSampler::sLogKernel(arma::uword b,
  arma::vec S_b,
  arma::vec cov_comb_log_det,
  arma::cube cov_comb_inv
) {
  
  double score = 0.0;
  arma::vec dist_from_mean(P);
  dist_from_mean.zeros();

  score = batchLikelihood(batch_ind(b), 
    labels, 
    cov_comb_log_det,
    t_df,
    mean_sum.cols(KB_inds + b),
    cov_comb_inv
  );
  
  for(arma::uword p = 0; p < P; p++) {
    score += -((rho + 1) * std::log(S_b(p) - S_loc) + theta / (S_b(p) - S_loc));
  }
  return score;
};

double mvtSampler::muLogKernel(arma::uword k, arma::vec mu_k, arma::mat mean_sum) {
  
  double score = 0.0;
  arma::uvec cluster_ind = arma::find(labels == k);
  arma::vec dist_from_mean(P);
  
  score = clusterLikelihood(
    t_df(k),
    cluster_ind,
    cov_comb_log_det.row(k).t(),
    mean_sum,
    cov_comb_inv.slices(k * B + B_inds)
  );
  score += -0.5 * arma::as_scalar(kappa * ((mu_k - xi).t() *  cov_inv.slice(k) * (mu_k - xi)));
  return score;
};

double mvtSampler::covLogKernel(arma::uword k, 
  arma::mat cov_k,
  double cov_log_det,
  arma::mat cov_inv,
  arma::vec cov_comb_log_det,
  arma::cube cov_comb_inv
) {
  
  double score = 0.0;
  arma::uvec cluster_ind = arma::find(labels == k);
  arma::vec dist_from_mean(P);
  
  score = clusterLikelihood(
    t_df(k),
    cluster_ind,
    cov_comb_log_det,
    mean_sum.cols(k * B + B_inds),
    cov_comb_inv
  );
  
  score += -0.5 *( arma::as_scalar((nu + P + 2) * cov_log_det 
      + kappa * ((mu.col(k) - xi).t() * cov_inv * (mu.col(k) - xi)) 
      + arma::trace(scale * cov_inv)
    )
  );

  return score;
};

double mvtSampler::dfLogKernel(arma::uword k, 
  double t_df,
  double pdf_coef
) {
  
  arma::uword b = 0;
  double score = 0.0;
  arma::uvec cluster_ind = arma::find(labels == k);
  arma::vec dist_from_mean(P);
  for (auto& n : cluster_ind) {
    b = batch_vec(n);
    dist_from_mean = X_t.col(n) - mean_sum.col(k * B + b);
    score += pdf_coef - 0.5 * (t_df + P) * log(1 + (1/t_df) * arma::as_scalar(dist_from_mean.t() * cov_comb_inv.slice(k * B + b) * dist_from_mean));
  }
  score += (psi - 1) * log(t_df - t_loc) - chi * (t_df - t_loc);
  return score;
};

void mvtSampler::clusterDFMetropolis() {
  double u = 0.0,
    proposed_model_score = 0.0, 
    acceptance_prob = 0.0, 
    current_model_score = 0.0, 
    t_df_proposed = 0.0, 
    proposed_pdf_coef = 0.0;
  
  for(arma::uword k = 0; k < K ; k++) {
    proposed_model_score = 0.0, acceptance_prob = 0.0, current_model_score = 0.0, t_df_proposed = 0.0, proposed_pdf_coef = 0.0;
    if(N_k(k) == 0){
      t_df_proposed = t_loc + arma::randg<double>( arma::distr_param(psi, 1.0 / chi));
      proposed_pdf_coef = calcPDFCoef(t_df_proposed);
    } else {
      
      // Propose a new value for the degrees of freedom
      t_df_proposed = t_loc + arma::randg( arma::distr_param( (t_df(k) - t_loc) * t_df_proposal_window, 1.0 / t_df_proposal_window) );
      
      // Proposed value
      proposed_pdf_coef = calcPDFCoef(t_df_proposed);
      
      // Asymmetric proposal density
      proposed_model_score = gammaLogLikelihood(t_df(k) - t_loc, (t_df_proposed - t_loc) * t_df_proposal_window, t_df_proposal_window);
      current_model_score = gammaLogLikelihood(t_df_proposed - t_loc, (t_df(k) - t_loc) * t_df_proposal_window, t_df_proposal_window);
      
      // The prior is included in the kernel
      proposed_model_score += dfLogKernel(k, t_df_proposed, proposed_pdf_coef);
      current_model_score += dfLogKernel(k, t_df(k), pdf_coef(k));
      
      u = arma::randu();
      acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));
      
    }
    
    if((u < acceptance_prob) || (N_k(k) == 0)) {
      t_df(k) = t_df_proposed;
      t_df_count(k)++;
      pdf_coef(k) = proposed_pdf_coef;
    }
  }
}

// Full-dataset Student-t log-likelihood (t_df held fixed; pdf_coef and
// cov_comb_log_det are constant w.r.t. gamma and dropped, as elsewhere in
// this class) under a given mean_sum, used by interactionMetropolis().
double mvtSampler::interactionDataLogLikelihood(arma::mat mean_sum_arg) {
  double score = 0.0;
  for(arma::uword n = 0; n < N; n++) {
    arma::uword c = labels(n) * B + batch_vec(n);
    arma::vec diff = X_t.col(n) - mean_sum_arg.col(c);
    double u = arma::as_scalar(diff.t() * cov_comb_inv.slice(c) * diff);
    score += -0.5 * (t_df(labels(n)) + P) * log(1.0 + (1.0 / t_df(labels(n))) * u);
  }
  return score;
};

// Block Metropolis-Hastings update for gamma(p, ., .) - see
// mvnSampler::interactionMetropolis() for the sum-to-zero-subspace
// derivation; the only difference here is the Student-t (rather than
// Gaussian) data log-likelihood used to score each candidate.
void mvtSampler::interactionMetropolis() {

  for(arma::uword p = 0; p < P; p++) {

    arma::mat raw_perturbation(K, B);
    for(arma::uword k = 0; k < K; k++) {
      for(arma::uword b = 0; b < B; b++) {
        raw_perturbation(k, b) = arma::randn() * gamma_proposal_window;
      }
    }
    arma::mat perturbation = doubleCenterMatrix(raw_perturbation);

    arma::mat proposed_mean_sum = mean_sum;
    double prior_current = 0.0, prior_proposed = 0.0;
    arma::mat current_gamma_p(K, B), proposed_gamma_p(K, B);

    for(arma::uword k = 0; k < K; k++) {
      for(arma::uword b = 0; b < B; b++) {
        current_gamma_p(k, b) = gamma(p, k, b);
        proposed_gamma_p(k, b) = current_gamma_p(k, b) + perturbation(k, b);

        prior_current += -0.5 * std::pow(current_gamma_p(k, b), 2.0) / tau2_interaction(p);
        prior_proposed += -0.5 * std::pow(proposed_gamma_p(k, b), 2.0) / tau2_interaction(p);

        proposed_mean_sum(p, k * B + b) = mean_sum(p, k * B + b) - current_gamma_p(k, b) + proposed_gamma_p(k, b);
      }
    }

    double current_model_score = interactionDataLogLikelihood(mean_sum) + prior_current;
    double proposed_model_score = interactionDataLogLikelihood(proposed_mean_sum) + prior_proposed;

    double u = arma::randu();
    double acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));

    if(u < acceptance_prob) {
      for(arma::uword k = 0; k < K; k++) {
        for(arma::uword b = 0; b < B; b++) {
          gamma(p, k, b) = proposed_gamma_p(k, b);
        }
      }
      mean_sum.row(p) = proposed_mean_sum.row(p);
      gamma_count(p)++;
    }
  }
};

void mvtSampler::metropolisStep() {

  // Metropolis step for cluster parameters
  clusterCovarianceMetropolis();
  clusterMeanMetropolis();
  clusterDFMetropolis();

  // Metropolis step for batch parameters
  if(sample_m_scale) {
    sampleMScalePosterior();
  }
  batchScaleMetropolis();
  batchShiftMetorpolis();

  if(include_interaction) {
    interactionMetropolis();
    sampleTauInteractionPosterior();
  }
};

