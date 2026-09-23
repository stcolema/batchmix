// mvnSamplerSeparationStrategy.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "pdfs.h"
# include "sampler.h"
# include "mvnSamplerSeparationStrategy.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// mvnSamplerSeparationStrategy class

mvnSamplerSeparationStrategy::mvnSamplerSeparationStrategy(
  arma::uword _K,
  arma::uword _B,
  double _mu_proposal_window,
  double _r_proposal_window,
  double _sigma_proposal_window,
  double _m_proposal_window,
  double _S_proposal_window,
  arma::uvec _labels,
  arma::uvec _batch_vec,
  arma::vec _concentration,
  arma::mat _X,
  arma::uvec _fixed,
  double _m_scale,
  double _rho,
  double _theta,
  bool _sample_m_scale,
  double _eta
) : sampler(_K,
_B,
_labels,
_batch_vec,
_concentration,
_X,
_fixed) {

  eta = _eta;

  // Empirical-Bayes prior hyperparameters below need a finite matrix -
  // mean()/cov() do not skip NaN - so use a column-mean-imputed copy for
  // this one-off calculation only. It has no bearing on the actual
  // per-sweep missing-data augmentation carried out by updateLatentData(),
  // which works from the untouched X_raw/items_to_augment (see sampler.h)
  // instead; X/X_t themselves are also given the same starting fill-in
  // below, purely as a safe placeholder until updateLatentData() first
  // runs (see the driver loop in sampleSemisupervisedMVNSeparationStrategy.cpp).
  mat X_imputed = imputeColumnMeans(X);

  rowvec X_min = min(X_imputed), X_max = max(X_imputed);
  mat global_cov = arma::cov(X_imputed);

  // Mean
  mat mean_mat = mean(X_imputed, 0).t();
  mu_0 = mean_mat.col(0);

  // Empirical Bayes fora diagonal covariance matrix
  mat scale_param = X_imputed.each_row() - mu_0.t();
  vec diag_entries(P);
  // double scale_entry = accu(scale_param % scale_param, 0) / (N * std::pow(K, 1.0 / (double) P));
  
  double scale_entry = (accu(global_cov.diag()) / P) / std::pow(K, 2.0 / (double) P);
  
  diag_entries.fill(scale_entry);
  scale = diagmat( diag_entries );
  
  scale = global_cov / std::pow(K, 2.0 / (double) P);
  
  // The mean of the prior distribution for the batch shift, m, parameter
  delta_2 = 0.0;
  lambda_2 = _m_scale;
  m_scale = _m_scale;
  
  // Prior precision is the inverse of something on the scale of 1/10 the global 
  // covariance
  delta_2 = accu(global_cov.diag()) / (double) P ;
  batch_shift_prior_precision = 1.0 / (delta_2 * lambda_2);
  
  t = 1.0 / ((accu(global_cov.diag()) / P ) * m_scale);
  
  // Hyperparameters for the batch scale
  rho = _rho;
  theta = _theta;

  // Set the size of the objects to hold the component specific parameters
  mu.set_size(P, K);
  mu.zeros();
  
  Rcpp::Rcout << "\nSetting new parameters sizes.\n";
  
  sigma.set_size(P, K);
  sigma.zeros();
  
  Sigma_mat.set_size(P, P, K);
  Sigma_mat.zeros();
  
  R.set_size(P, P, K);
  R.zeros();
  
  cov.set_size(P, P, K);
  cov.zeros();
  
  // Set the size of the objects to hold the batch specific parameters
  m.set_size(P, B);
  m.zeros();
  
  // We are assuming a diagonal structure in the batch scale
  S.set_size(P, B);
  S.zeros();
  
  // Count the number of times proposed values are accepted
  r_count = zeros<uvec>(K);
  sigma_count = zeros<uvec>(K);
  mu_count = zeros<uvec>(K);
  m_count = zeros<uvec>(B);
  S_count = zeros<uvec>(B);
  
  // These will hold vertain matrix operations to avoid computational burden
  // The log determinant of each cluster covariance
  r_log_det = zeros<vec>(K);
  cov_log_det = zeros<vec>(K);
  
  // The log determinant of the covariance combination
  cov_comb_log_det.set_size(K, B);
  cov_comb_log_det.zeros();
  
  // The possible combinations for the sum of the cluster and batch means
  mean_sum.set_size(P, K * B);
  mean_sum.zeros();
  
  // The combination of each possible cluster and batch covariance
  cov_comb.set_size(P, P, K * B);
  cov_comb.zeros();
  
  // Inverse of the cluster covariance
  cov_inv.set_size(P, P, K);
  cov_inv.zeros();
  
  // The inverse of the covariance combination
  cov_comb_inv.set_size(P, P, K * B);
  cov_comb_inv.zeros();
  
  Rcpp::Rcout << "\nProposal windows.\n";
  
  // The batch corrected data
  Y.set_size(N, P);
  Y.zeros();

  // Working complete-data copy: seed missing (NaN) entries of X/X_t with
  // the observed column mean as a safe placeholder. updateLatentData()
  // (called at the top of every MCMC sweep - see the driver loop in
  // sampleSemisupervisedMVNSeparationStrategy.cpp) replaces this with a
  // proper draw from the full conditional before it is ever used in a
  // likelihood; this fill-in only guards against reading uninitialised/NaN
  // memory before that first call.
  for(uword p = 0; p < P; p++) {
    for(uword n = 0; n < N; n++) {
      if(!std::isfinite(X(n, p))) {
        X(n, p) = X_imputed(n, p);
      }
    }
  }
  X_t = X.t();

  // The proposal windows for the cluster and batch parameters
  mu_proposal_window = _mu_proposal_window;
  r_proposal_window = _r_proposal_window;
  sigma_proposal_window = _sigma_proposal_window;
  m_proposal_window = _m_proposal_window;
  S_proposal_window = _S_proposal_window;
  
  // flag indicating if the hyperparameter on the batch shift is sampled or not
  sample_m_scale = _sample_m_scale;
};

void mvnSamplerSeparationStrategy::sampleCovPrior() {
  vec log_sigma(P);
  log_sigma.zeros();
  mat cov_check;
  for(uword k = 0; k < K; k++){
    // R/sigma are each individually valid by construction (LKJ always
    // gives a PD correlation matrix; exp() of a normal draw is always
    // positive), but their product can, with vanishingly small
    // probability, be PD in exact arithmetic yet numerically singular in
    // floating point once R sits close to the boundary of the PD cone
    // (confirmed empirically - see sigmaMHStep()/rMHStep()'s equivalent
    // comments). matrixCombinations() runs immediately after
    // sampleFromPriors() and unconditionally inverts cov.slice(k), so an
    // initial draw this pathological would crash the sampler before a
    // single sweep ran; redraw (a negligible-measure rejection of the
    // pathological tail, not a bias on the prior) rather than let that
    // happen.
    bool ok = false;
    const uword max_attempts = 50;
    for(uword attempt = 0; attempt < max_attempts && !ok; attempt++) {
      R.slice(k) = sampleLKJCorrelationMatrix(P, eta);
      log_sigma = randn<vec>(P, distr_param(beta, xi));

      sigma.col(k) = exp(log_sigma);
      Sigma_mat.slice(k).diag() = sigma.col(k);

      cov.slice(k) = Sigma_mat.slice(k) * R.slice(k) * Sigma_mat.slice(k);
      ok = arma::inv_sympd(cov_check, cov.slice(k));
    }
  }
};

void mvnSamplerSeparationStrategy::sampleMuPrior() {
  for(uword k = 0; k < K; k++){
    mu.col(k) = mvnrnd(mu_0, (1.0/kappa) * cov.slice(k), 1);
  }
};

void mvnSamplerSeparationStrategy::sampleSPrior() {
  for(uword b = 0; b < B; b++){
    for(uword p = 0; p < P; p++){
      S(p, b) = S_loc + 1.0 / randg<double>( distr_param(rho, 1.0 / theta ) );
    }
  }
};

void mvnSamplerSeparationStrategy::sampleMPrior() {
  // batch_shift_prior_precision is a precision (1 / (delta_2 * lambda_2)),
  // so the prior standard deviation is its inverse square root, not the
  // precision itself.
  double batch_shift_prior_sd = std::sqrt(1.0 / batch_shift_prior_precision);
  for(uword b = 0; b < B; b++){
    for(uword p = 0; p < P; p++){
      m(p, b) = randn<double>() * batch_shift_prior_sd + batch_shift_prior_mean;
    }
  }
};

void mvnSamplerSeparationStrategy::sampleFromPriors() {
  if(sample_m_scale) {
    sampleMScalePrior();
  }
  sampleCovPrior();
  sampleMuPrior();
  sampleSPrior();
  sampleMPrior();
  if(include_interaction) {
    sampleGammaPrior();
  }
};

void mvnSamplerSeparationStrategy::sampleMScalePrior() {
  lambda_2 = rInvGamma(a, b);
  batch_shift_prior_precision = 1.0 / (delta_2 * lambda_2);
}

void mvnSamplerSeparationStrategy::sampleMScalePosterior() {
  double a_pos = 0.0, b_pos = 0.0;
  a_pos = a + 0.5 * P * B;
  b_pos = accu(pow(m, 2.0)) / (2.0 * delta_2) + b;
  lambda_2 = rInvGamma(a_pos, b_pos);
  batch_shift_prior_precision = 1.0 / (delta_2 * lambda_2);
}

// Update the common matrix manipulations to avoid recalculating N times
void mvnSamplerSeparationStrategy::matrixCombinations() {
  
  for(uword k = 0; k < K; k++) {
    cov_inv.slice(k) = inv_sympd(cov.slice(k));
    cov_log_det(k) = log_det(cov.slice(k)).real();
    for(uword b = 0; b < B; b++) {
      cov_comb.slice(k * B + b) = cov.slice(k);
      for(uword p = 0; p < P; p++) {
        cov_comb.slice(k * B + b)(p, p) *= S(p, b);
      }
      cov_comb_log_det(k, b) = log_det(cov_comb.slice(k * B + b)).real();
      cov_comb_inv.slice(k * B + b) = inv_sympd(cov_comb.slice(k * B + b));

      mean_sum.col(k * B + b) = mu.col(k) + m.col(b);
      if(include_interaction) {
        mean_sum.col(k * B + b) += gamma.slice(b).col(k);
      }
    }
  }
};

// See the header comment for the provenance/derivation; this is a direct
// port of mvnSamplerMixed::updateLatentData() with the binary/censoring
// branches removed.
void mvnSamplerSeparationStrategy::updateLatentData() {

  vec z_i(P), eta_mean(P);
  double cond_mean = 0.0, cond_var = 0.0, cond_sd = 0.0, lambda_pp = 0.0;
  uword k = 0, b = 0, kb = 0;

  for(auto& n : items_to_augment) {

    k = labels(n);
    b = batch_vec(n);
    kb = k * B + b;

    z_i = X_t.col(n);
    eta_mean = mean_sum.col(kb);

    for(uword p = 0; p < P; p++) {

      if(std::isfinite(X_raw_t(p, n))) {
        continue;
      }

      lambda_pp = cov_comb_inv(p, p, kb);
      cond_var = 1.0 / lambda_pp;
      cond_mean = eta_mean(p) - cond_var * (
        dot(cov_comb_inv.slice(kb).row(p), z_i - eta_mean) - lambda_pp * (z_i(p) - eta_mean(p))
      );
      cond_sd = std::sqrt(cond_var);

      z_i(p) = cond_mean + cond_sd * randn();
    }

    X_t.col(n) = z_i;
  }

  X = X_t.t();
};

// The log likelihood of a item belonging to each cluster given the batch label.
arma::vec mvnSamplerSeparationStrategy::itemLogLikelihood(arma::vec x, arma::uword b) {
  
  double exponent = 0.0;
  vec ll(K), dist_to_mean(P), m_b(B);
  ll.zeros();
  dist_to_mean.zeros();
  m_b = m.col(b);
  
  for(uword k = 0; k < K; k++){
    
    // The exponent part of the MVN pdf
    dist_to_mean = x - mean_sum.col(k * B + b);
    exponent = as_scalar(dist_to_mean.t() * cov_comb_inv.slice(k * B + b) * dist_to_mean);
    
    // Normal log likelihood
    ll(k) = -0.5 *(cov_comb_log_det(k, b) + exponent + (double) P * log(2.0 * M_PI));
  }
  
  return(ll);
};

void mvnSamplerSeparationStrategy::calcBIC(){

  // Each occupied component has a weight, a mean vector and a symmetric
  // covariance matrix (via its R/sigma decomposition); each batch has a
  // shift vector and a scale vector. structuralExtraBICParams() adds the
  // interaction term's and/or the batch-specific weight prior's extra
  // parameters when either is enabled (0 in the default configuration).
  BIC = 2 * observed_likelihood - (n_param_cluster * K_occ + n_param_batch * B + structuralExtraBICParams()) * std::log(N);

};

double mvnSamplerSeparationStrategy::groupLikelihood(arma::uvec inds,
                               arma::uvec group_inds,
                               arma::vec cov_det,
                               arma::mat mean_sum,
                               arma::cube cov_inv){
  
  // Rcpp::Rcout << "\nCov det:\n" << cov_det.t();
  // Rcpp::Rcout << "\n\nMean sum:\n" << mean_sum;
  // Rcpp::Rcout << "\n\nCov inverse:\n" << cov_inv;
  // Rcpp::Rcout << "\n\nBatches:\n" << unique(group_inds).t();
  
  uword c = 0;
  double score = 0.0;
  vec dist_from_mean(P);
  
  for (auto& n : inds) {
    c = group_inds(n);
    dist_from_mean = X_t.col(n) - mean_sum.col(c);
    score += as_scalar(cov_det(c) + (dist_from_mean.t() * cov_inv.slice(c) * dist_from_mean));
  }
  return (-0.5 * score);
}


double mvnSamplerSeparationStrategy::mLogKernel(arma::uword b, arma::vec m_b, arma::mat mean_sum) {
  
  double score = 0.0;
  vec dist_from_mean(P);
  dist_from_mean.zeros();
  
  score = groupLikelihood(batch_ind(b),
    labels,
    cov_comb_log_det.col(b),
    mean_sum,
    cov_comb_inv.slices(KB_inds + b)
  );
  
  for(uword p = 0; p < P; p++) {
    score += -0.5 * (batch_shift_prior_precision * std::pow(m_b(p) - batch_shift_prior_mean, 2.0) );
  }

  return score;
};

double mvnSamplerSeparationStrategy::sLogKernel(arma::uword b, 
                          arma::vec S_b, 
                          arma::vec cov_comb_log_det,
                          arma::cube cov_comb_inv) {
  
  double score = 0.0;
  vec dist_from_mean(P);
  dist_from_mean.zeros();
  
  score = groupLikelihood(batch_ind(b),
    labels,
    cov_comb_log_det,
    mean_sum.cols(KB_inds + b),
    cov_comb_inv
  );
  
  for(uword p = 0; p < P; p++) {
    score +=  -((rho + 1) * std::log(S_b(p) - S_loc) + theta / (S_b(p) - S_loc));
  }
  return score;
};

double mvnSamplerSeparationStrategy::muLogKernel(arma::uword k, arma::vec mu_k, arma::mat mean_sum) {
  
  double score = 0.0;
  uvec cluster_ind = find(labels == k);
  vec dist_from_mean(P);
  
  score = groupLikelihood(cluster_ind,
    batch_vec,
    cov_comb_log_det.row(k).t(),
    mean_sum,
    cov_comb_inv.slices(k * B + B_inds)
  );
  
  score += -0.5 * as_scalar(kappa * ((mu_k - mu_0).t() *  cov_inv.slice(k) * (mu_k - mu_0)));
  
  return score;
};


double mvnSamplerSeparationStrategy::rLogKernel(uword k,
                              double r_log_det,
                              vec cov_comb_log_det,
                              mat cov,
                              mat cov_inverse,
                              cube cov_comb_inv
)
{
  double score = 0.0, exponent = 0.0;
  uvec cluster_ind = find(labels == k);
  vec dist_from_prior = mu.col(k) - mu_0;

  score = groupLikelihood(cluster_ind,
    batch_vec,
    cov_comb_log_det,
    mean_sum.cols(k * B + B_inds),
    cov_comb_inv
  );

  exponent = arma::as_scalar(dist_from_prior.t() * cov_inverse * dist_from_prior);

  // LKJ(eta) prior on R (unnormalised; the normalising constant depends
  // only on eta and P, not on R, so it cancels in the Metropolis ratio),
  // plus the contribution of the mu | cov conditional prior that depends
  // on cov through R.
  score += (eta - 1.0) * r_log_det
    -0.5 * kappa * exponent
    -0.5 * log_det((1.0 / kappa) * cov).real();

  return score;
};

double mvnSamplerSeparationStrategy::sigmaLogKernel(uword k,
                                  vec cov_comb_log_det,
                                  mat sigma_mat,
                                  mat cov,
                                  mat cov_inverse,
                                  cube cov_comb_inv
) {
  double score = 0.0, exponent = 0.0;
  uvec cluster_ind = find(labels == k);
  arma::vec dist_from_prior = mu.col(k) - mu_0;
  
  score = groupLikelihood(cluster_ind,
                          batch_vec,
                          cov_comb_log_det,
                          mean_sum.cols(k * B + B_inds),
                          cov_comb_inv
  );
  exponent = arma::as_scalar(dist_from_prior.t() * cov_inverse * dist_from_prior);

  // Log-normal prior: log(sigma_{k,p}) ~ Normal(beta, xi^2) (matching how
  // it is drawn in sampleCovPrior: randn(P, distr_param(beta, xi))).
  // Evaluated at the candidate/current value actually passed in
  // (sigma_mat), NOT the class member sigma (which, critically, is the
  // *other* cluster's value unless k == 0, and in any case is not the
  // value being scored by this call - using it here was the bug).
  for(uword p = 0; p < P; p++) {
    score += -std::pow(log(sigma_mat(p)) - beta, 2.0) / (2.0 * xi * xi);
  }
  
  score += -0.5 * log_det((1.0 / kappa) * cov).real()
    -0.5 * kappa * exponent;
  
  // Rcpp::Rcout << "SIGMA KERNEL: Finished.\n\n";
  
  return score;
};

void mvnSamplerSeparationStrategy::batchScaleMetropolis() {
  
  bool next = false;
  double u = 0.0, proposed_model_score = 0.0, acceptance_prob = 0.0, current_model_score = 0.0;
  vec S_proposed(P), proposed_cov_comb_log_det(K);
  cube proposed_cov_comb(P, P, K), proposed_cov_comb_inv(P, P, K);
  
  S_proposed.zeros();
  proposed_cov_comb_log_det.zeros();
  proposed_cov_comb.zeros();
  proposed_cov_comb_inv.zeros();
  
  for(uword b = 0; b < B; b++) {
    
    next = false;
    acceptance_prob = 0.0, proposed_model_score = 0.0, current_model_score = 0.0;
    proposed_cov_comb.zeros();
    
    for(uword p = 0; p < P; p++) {
      if((S(p, b) - S_loc) * S_proposal_window < 0.0){ 
        Rcpp::stop("\n\nCurent batch scale equals S_loc");
      }
      if( (1.0 / S_proposal_window) < 0.0) {
        Rcpp::stop("\n\nCurent batch scale proposal window is some how negative?");
      }
      S_proposed(p) = S_loc + randg( distr_param( (S(p, b) - S_loc) * S_proposal_window, 1.0 / S_proposal_window) );
      
      if(S_proposed(p) <= S_loc) {
        next = true;
      }

      // Metropolis-Hastings correction for the asymmetric Gamma proposal.
      // proposed_model_score must accumulate log pi(proposed) + log q(current | proposed)
      // (the REVERSE proposal density) and current_model_score must
      // accumulate log pi(current) + log q(proposed | current) (the FORWARD
      // density, i.e. the density of the draw actually taken). This was
      // previously swapped here (and in every other copy of this Gamma-RW
      // pattern in the package bar mvtSampler::clusterDFMetropolis), which
      // silently biases the chain - verified by simulating this exact
      // proposal against a known Gamma(5, 2) target, where the swapped
      // assignment recovers a posterior mean of ~1.41 instead of the true 2.5.
      current_model_score += gammaLogLikelihood(S_proposed(p) - S_loc, (S(p, b) - S_loc) * S_proposal_window, S_proposal_window);
      proposed_model_score += gammaLogLikelihood(S(p, b) - S_loc, (S_proposed(p) - S_loc) * S_proposal_window, S_proposal_window);
    }
    
    if(next) {
      continue;
    }
    
    proposed_cov_comb = cov;
    // S_proposed > S_loc always (a Gamma draw shifted by S_loc), so this
    // only inflates each cluster's diagonal, which alone keeps PD +
    // PSD = PD exactly - but cov.slice(k) itself may already be close to
    // the boundary of the PD cone (see rMHStep()'s/sigmaMHStep()'s
    // equivalent comments), so guard here too rather than assume.
    // Automatic-reject-on-degenerate-proposal, same as those.
    bool cov_ok = true;
    for(uword k = 0; k < K; k++) {
      for(uword p = 0; p < P; p++) {
        proposed_cov_comb.slice(k)(p, p) *= S_proposed(p);
      }
      proposed_cov_comb_log_det(k) = log_det(proposed_cov_comb.slice(k)).real();
      cov_ok = arma::inv_sympd(proposed_cov_comb_inv.slice(k), proposed_cov_comb.slice(k));
      if(!cov_ok) {
        break;
      }
    }
    if(!cov_ok) {
      continue;
    }

    proposed_model_score += sLogKernel(b,
      S_proposed, 
      proposed_cov_comb_log_det,
      proposed_cov_comb_inv
    );
    
    current_model_score += sLogKernel(b, 
      S.col(b), 
      cov_comb_log_det.col(b),
      cov_comb_inv.slices(KB_inds + b)
    );
    
    u = randu();
    acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));
    
    if(u < acceptance_prob){
      S.col(b) = S_proposed;
      S_count(b)++;
      
      for(uword k = 0; k < K; k++) {
        cov_comb.slice(k * B + b) = proposed_cov_comb.slice(k);
        cov_comb_log_det(k, b) = proposed_cov_comb_log_det(k);
        cov_comb_inv.slice(k * B + b) = proposed_cov_comb_inv.slice(k);
      }
    }
  }
};

void mvnSamplerSeparationStrategy::batchShiftMetorpolis() {
  
  double u = 0.0, 
    proposed_model_score = 0.0, 
    acceptance_prob = 0.0, 
    current_model_score = 0.0;
  
  vec m_proposed(P);
  mat proposed_mean_sum(P, K);
  m_proposed.zeros();
  
  for(arma::uword b = 0; b < B; b++) {
    for(arma::uword p = 0; p < P; p++){
      // The proposal window is now a diagonal matrix of common entries.
      m_proposed(p) = (arma::randn() * m_proposal_window) + m(p, b);
    }
    
    for(arma::uword k = 0; k < K; k++) {
      proposed_mean_sum.col(k) = mu.col(k) + m_proposed;
      if(include_interaction) {
        proposed_mean_sum.col(k) += gamma.slice(b).col(k);
      }
    }

    // We have a symmetric proposal density so the posterior kernel is the only
    // thing we're interested in
    proposed_model_score = mLogKernel(b, m_proposed, proposed_mean_sum);
    current_model_score = mLogKernel(b, m.col(b), mean_sum.cols(KB_inds + b));
    

    u = arma::randu();
    acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));
    
    // If we accept the proposed value, update our count of acceptances, the 
    // batch shift and the mean sum
    if(u < acceptance_prob){
      m.col(b) = m_proposed;
      m_count(b)++;
      
      for(arma::uword k = 0; k < K; k++) {
        mean_sum.col(k * B + b) = proposed_mean_sum.col(k);
      }
    }
  }
};

void mvnSamplerSeparationStrategy::rMHStep() {
  double u = 0.0,
    proposed_model_score = 0.0,
    acceptance_prob = 0.0,
    current_model_score = 0.0,
    proposed_r_log_det = 0.0,
    proposed_cov_log_det = 0.0,
    log_jacobian_ratio = 0.0;

  arma::vec proposed_cov_comb_log_det(B);
  arma::mat R_proposed(P, P),
    proposed_cov(P, P),
    proposed_cov_inv(P, P);
  arma::cube proposed_cov_comb(P, P, B), proposed_cov_comb_inv(P, P, B);

  R_proposed.zeros();
  proposed_cov.zeros();
  proposed_cov_inv.zeros();
  proposed_cov_comb_log_det.zeros();
  proposed_cov_comb.zeros();
  proposed_cov_comb_inv.zeros();

  for(arma::uword k = 0; k < K ; k++) {

    proposed_cov_comb.zeros();

    acceptance_prob = 0.0, proposed_model_score = 0.0, current_model_score = 0.0;

    // If no items in the class, sample from the prior distribution.
    //
    // R_proposed is a fresh, unconstrained LKJ(eta) draw, combined here
    // with Sigma_mat.slice(k) - the CURRENT sigma, which sigmaMHStep()'s
    // own empty-cluster branch may itself have freshly redrawn from its
    // LogNormal(beta, xi) prior on a previous sweep (see the comment
    // there). Neither draw is extreme on its own, but the combination is,
    // with vanishingly small but nonzero probability over many thousands
    // of sweeps, PD in exact arithmetic yet numerically singular in
    // floating point (confirmed empirically) - inv_sympd() throwing on
    // that is a real, reproducible crash. Redraw R_proposed (a
    // negligible-measure rejection of the pathological tail, not a bias
    // on the prior) up to max_attempts times using the non-throwing
    // two-argument inv_sympd(), and if every attempt is still numerically
    // degenerate, leave this cluster's R/cov untouched for this sweep
    // rather than crash.
    if(N_k(k) == 0){
      bool ok = false;
      const arma::uword max_attempts = 50;
      for(arma::uword attempt = 0; attempt < max_attempts && !ok; attempt++) {

        R_proposed = sampleLKJCorrelationMatrix(P, eta);
        proposed_cov = Sigma_mat.slice(k) * R_proposed * Sigma_mat.slice(k);
        ok = arma::inv_sympd(proposed_cov_inv, proposed_cov);
        if(!ok) {
          continue;
        }

        for(arma::uword b = 0; b < B; b++) {
          proposed_cov_comb.slice(b) = proposed_cov;
          for(arma::uword p = 0; p < P; p++) {
            proposed_cov_comb.slice(b)(p, p) *= S(p, b);
          }
          ok = arma::inv_sympd(proposed_cov_comb_inv.slice(b), proposed_cov_comb.slice(b));
          if(!ok) {
            break;
          }
        }
      }

      if(!ok) {
        continue;
      }

      proposed_r_log_det = arma::log_det(R_proposed).real();
      proposed_cov_log_det = arma::log_det(proposed_cov).real();
      for(arma::uword b = 0; b < B; b++) {
        proposed_cov_comb_log_det(b) = arma::log_det(proposed_cov_comb.slice(b)).real();
      }
    } else {

      // Random-walk proposal via the Cholesky-factor/partial-correlation
      // reparameterisation (buildCorrelationCholeskyFromZ /
      // choleskyToPartialCorrelations / logJacobianZToR in
      // genericFunctions.h), rather than a direct random walk on R's
      // pairwise entries. Any candidate Z in (-1,1)^(P(P-1)/2) yields a
      // valid PD correlation matrix, so no positive-definiteness
      // rejection is ever needed here (unlike the previous pairwise
      // Fisher-z approach) - and the geometry is far better behaved:
      // random-walk proposals directly on pairwise correlations do not
      // respect the curved, bounded geometry of the correlation-matrix
      // manifold, and were found (via extensive diagnosis - see the
      // package's development notes) to leave the *joint* R/sigma chain
      // stuck in a self-consistent but biased region indefinitely, even
      // though each 1-D conditional was independently verified exactly
      // correct. This is the standard fix (Stan Reference Manual,
      // "Cholesky Factors of Correlation Matrices"; McElreath,
      // *Statistical Rethinking*, ch. 14).
      // arma::chol()'s 2-argument (output-parameter) form is the
      // non-throwing one: it returns false, rather than throwing, on
      // failure, and on failure leaves L_current EMPTY (0x0) rather than
      // P x P. The return value was not being checked here, so a failure
      // silently fed an empty matrix into choleskyToPartialCorrelations(),
      // which indexes it assuming P x P - an out-of-bounds crash
      // (confirmed empirically via gdb backtrace), not a hypothetical.
      // R.slice(k) is always a valid correlation matrix in exact
      // arithmetic (every acceptance in this file now checks inv_sympd()
      // succeeds first - see the comments above), but "invertible enough
      // for inv_sympd" and "well-conditioned enough for a stable Cholesky
      // factorisation" are not quite the same numerical threshold, so a
      // just-barely-accepted R can still fail here. Treat that the same
      // way every other numerically-degenerate case in this file is
      // treated: skip this cluster's correlation update for this sweep
      // (no candidate proposal can be formed, so there is nothing to
      // accept or reject) rather than crash.
      mat L_current;
      bool chol_ok = arma::chol(L_current, R.slice(k), "lower");
      if(!chol_ok) {
        continue;
      }
      mat Z_current = choleskyToPartialCorrelations(L_current, P);
      mat Z_proposed = Z_current;

      for(uword i = 1; i < P; i++) {
        for(uword j = 0; j < i; j++) {
          double w_ij_proposed = std::atanh(Z_current(i, j)) + arma::randn() * r_proposal_window;
          Z_proposed(i, j) = std::tanh(w_ij_proposed);
        }
      }

      mat L_proposed = buildCorrelationCholeskyFromZ(Z_proposed, P);
      R_proposed = L_proposed * L_proposed.t();

      // Full log-Jacobian: the composition of the tanh unconstraining
      // layer (per free partial correlation) and the partial-correlation-
      // to-R layer.
      double tanh_jac_current = 0.0, tanh_jac_proposed = 0.0;
      for(uword i = 1; i < P; i++) {
        for(uword j = 0; j < i; j++) {
          tanh_jac_current += std::log(1.0 - Z_current(i, j) * Z_current(i, j));
          tanh_jac_proposed += std::log(1.0 - Z_proposed(i, j) * Z_proposed(i, j));
        }
      }

      log_jacobian_ratio = (logJacobianZToR(Z_proposed, P) + tanh_jac_proposed)
        - (logJacobianZToR(Z_current, P) + tanh_jac_current);

      proposed_model_score = log_jacobian_ratio;

      proposed_cov = Sigma_mat.slice(k) * R_proposed * Sigma_mat.slice(k);
      // The Cholesky/partial-correlation reparameterisation guarantees
      // R_proposed is PD in exact arithmetic for any Z in (-1,1)^(P(P-1)/2)
      // - but it does not keep R_proposed away from the boundary of the PD
      // cone, and a proposal legitimately close to that boundary (e.g. two
      // dimensions with a highly correlated random-walk state) can be PD
      // in exact arithmetic yet numerically singular in floating point
      // (confirmed empirically: a real proposal here had eigenvalues down
      // to ~1.4e-4, with the smallest nominally -1.3e-16 after symmetrising
      // - a genuine near-singular case, not a coding bug in the
      // reparameterisation). The textbook-correct treatment of a proposal
      // whose target density is numerically undefined is to treat it as
      // automatically rejected (acceptance probability 0), exactly like
      // the `next`/`continue` pattern already used elsewhere in this file
      // for degenerate individual-parameter proposals - not to retry with
      // a different proposal, which is only appropriate for the
      // unconditional prior draw used in the N_k(k) == 0 branch above.
      bool cov_ok = arma::inv_sympd(proposed_cov_inv, proposed_cov);
      if(!cov_ok) {
        continue;
      }
      proposed_r_log_det = arma::log_det(R_proposed).real();
      proposed_cov_log_det = arma::log_det(proposed_cov).real();

      for(arma::uword b = 0; b < B; b++) {
        proposed_cov_comb.slice(b) = proposed_cov;
        for(arma::uword p = 0; p < P; p++) {
          proposed_cov_comb.slice(b)(p, p) *= S(p, b);
        }
        proposed_cov_comb_log_det(b) = arma::log_det(proposed_cov_comb.slice(b)).real();
        cov_ok = arma::inv_sympd(proposed_cov_comb_inv.slice(b), proposed_cov_comb.slice(b));
        if(!cov_ok) {
          break;
        }
      }
      if(!cov_ok) {
        continue;
      }

      // The proposed model score is updated by the posterior kernel with the
      // proposed covariance, it's log determinant, inverse and the new
      // likelihood covariance
      proposed_model_score += rLogKernel(k,
                                         proposed_r_log_det,
                                         proposed_cov_comb_log_det,
                                         proposed_cov,
                                         proposed_cov_inv,
                                         proposed_cov_comb_inv
      );

      // The current model score is based on the posterior kernel for the current
      // covariance matrix
      current_model_score = rLogKernel(k,
                                        r_log_det(k),
                                        cov_comb_log_det.row(k).t(),
                                        cov.slice(k),
                                        cov_inv.slice(k),
                                        cov_comb_inv.slices(k * B + B_inds)
      );

      // Accept or reject
      u = arma::randu();
      acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));

    }
    if( (u < acceptance_prob) || (N_k(k) == 0) ){
      r_count(k)++;
      
      R.slice(k) = R_proposed;
      r_log_det(k) = proposed_r_log_det;
      
      cov.slice(k) = proposed_cov;
      cov_inv.slice(k) = proposed_cov_inv;
      cov_log_det(k) = proposed_cov_log_det;
      for(arma::uword b = 0; b < B; b++) {
        cov_comb.slice(k * B + b) = proposed_cov_comb.slice(b);
        cov_comb_log_det(k, b) = proposed_cov_comb_log_det(b);
        cov_comb_inv.slice(k * B + b) = proposed_cov_comb_inv.slice(b);
      }
    }
  }
};

void mvnSamplerSeparationStrategy::sigmaMHStep() {
  bool next = false;
  
  double u = 0.0, 
    proposed_model_score = 0.0, 
    acceptance_prob = 0.0, 
    current_model_score = 0.0, 
    proposed_sigma_log_det = 0.0,
    proposed_cov_log_det = 0.0;
  
  arma::vec sigma_proposed(P), proposed_cov_comb_log_det(B);
  arma::mat sigma_mat_proposed(P, P),proposed_cov(P, P), proposed_cov_inv(P, P);
  arma::cube proposed_cov_comb(P, P, B), proposed_cov_comb_inv(P, P, B);
  
  sigma_proposed.zeros();
  sigma_mat_proposed.zeros();
  
  proposed_cov.zeros();
  proposed_cov_inv.zeros();
  proposed_cov_comb_log_det.zeros();
  proposed_cov_comb.zeros();
  proposed_cov_comb_inv.zeros();
  
  for(arma::uword k = 0; k < K ; k++) {
    
    // Rcpp::Rcout << "\n\nSIGMA: Loop " << k << "\n";
    
    sigma_mat_proposed.diag() = sigma.col(k);
    next = false;
    
    proposed_cov_comb.zeros();
    
    acceptance_prob = 0.0, proposed_model_score = 0.0, current_model_score = 0.0;
    
    // If no items in the class, sample fresh from the prior distribution
    // (matching sampleCovPrior()'s log_sigma ~ N(beta, xi), sigma = exp(log_sigma))
    // rather than an uncorrected, force-accepted random walk off the
    // cluster's current (possibly stale) sigma - see clusterMeanMetropolis()'s
    // N_k(k)==0 branch for the same pattern applied to mu.
    //
    // A fresh LogNormal(beta, xi) draw is, with vanishingly small but
    // nonzero probability over many thousands of sweeps, extreme enough
    // (xi = 1 on the log scale) that sigma_mat_proposed * R.slice(k) *
    // sigma_mat_proposed is PD in exact arithmetic but numerically
    // singular in floating point once R sits close to the boundary of the
    // PD cone (LKJ(eta) does not keep correlations away from that
    // boundary) - inv_sympd() throwing on that is a real, reproducible
    // crash (confirmed empirically), not a hypothetical. The single-arg
    // inv_sympd() used everywhere else in this file assumes that never
    // happens; here, where the whole point is an unconstrained fresh draw
    // every sweep an empty cluster persists, redraw (a negligible-measure
    // rejection of the pathological tail, not a bias on the prior) up to
    // max_attempts times using the non-throwing two-argument inv_sympd(),
    // and if every attempt is still numerically degenerate, leave this
    // cluster's sigma/cov untouched for this sweep rather than crash - the
    // same "reject this sweep's proposal" behaviour degenerate/singular
    // proposals already get via `next`/`continue` in the branch below.
    if(N_k(k) == 0){

      bool ok = false;
      const arma::uword max_attempts = 50;
      for(arma::uword attempt = 0; attempt < max_attempts && !ok; attempt++) {

        sigma_proposed = arma::exp(arma::randn<vec>(P, distr_param(beta, xi)));
        sigma_mat_proposed.diag() = sigma_proposed;

        proposed_cov = sigma_mat_proposed * R.slice(k) * sigma_mat_proposed;
        ok = arma::inv_sympd(proposed_cov_inv, proposed_cov);
        if(!ok) {
          continue;
        }

        for(arma::uword b = 0; b < B; b++) {
          proposed_cov_comb.slice(b) = proposed_cov;
          for(arma::uword p = 0; p < P; p++) {
            proposed_cov_comb.slice(b)(p, p) *= S(p, b);
          }
          ok = arma::inv_sympd(proposed_cov_comb_inv.slice(b), proposed_cov_comb.slice(b));
          if(!ok) {
            break;
          }
        }
      }

      if(!ok) {
        // Every attempt was numerically degenerate - leave this cluster's
        // parameters untouched this sweep rather than propagate a matrix
        // inv_sympd() could not invert.
        continue;
      }

      proposed_sigma_log_det = arma::log_det(sigma_mat_proposed).real();
      proposed_cov_log_det = arma::log_det(proposed_cov).real();
      for(arma::uword b = 0; b < B; b++) {
        proposed_cov_comb_log_det(b) = arma::log_det(proposed_cov_comb.slice(b)).real();
      }
    } else {

      // Proposal density is a Wishart with an expected value of the current
      // covariance matrix
      for(uword p = 0; p < P; p++) {

        // Rcpp::Rcout << "Inter Loop: " << p << "\n";

        sigma_proposed(p) = randg( distr_param( sigma(p, k) * sigma_proposal_window, 1.0 / sigma_proposal_window) );

        if(sigma_proposed(p) <= 1e-8) {
          next = true;
        }

        // Asymmetric proposal density (see the derivation/empirical check
        // above): reverse density q(current|proposed) -> proposed_model_score,
        // forward density q(proposed|current) -> current_model_score.
        proposed_model_score += gammaLogLikelihood(sigma(p, k), sigma_proposed(p) * sigma_proposal_window, sigma_proposal_window);
        current_model_score += gammaLogLikelihood(sigma_proposed(p), sigma(p, k) * sigma_proposal_window, sigma_proposal_window);
      }
      
      if(next) {
        continue;
      }
      
      // Rcpp::Rcout << "\n\nSIGMA MH: Value proposed.\n";
      
      sigma_mat_proposed.diag() = sigma_proposed;
      
      proposed_cov = sigma_mat_proposed * R.slice(k) * sigma_mat_proposed;
      // Same reasoning as rMHStep()'s non-empty branch: a proposal can be
      // PD in exact arithmetic (R.slice(k) is always a valid correlation
      // matrix, sigma_proposed always positive) yet numerically singular
      // in floating point when R.slice(k) sits close to the boundary of
      // the PD cone - treat that as an automatic reject (acceptance
      // probability 0), not a retry.
      bool cov_ok = arma::inv_sympd(proposed_cov_inv, proposed_cov);
      if(!cov_ok) {
        continue;
      }
      proposed_sigma_log_det = arma::log_det(sigma_mat_proposed).real();
      proposed_cov_log_det = arma::log_det(proposed_cov).real();

      for(arma::uword b = 0; b < B; b++) {
        proposed_cov_comb.slice(b) = proposed_cov;
        for(arma::uword p = 0; p < P; p++) {
          proposed_cov_comb.slice(b)(p, p) *= S(p, b);
        }
        proposed_cov_comb_log_det(b) = arma::log_det(proposed_cov_comb.slice(b)).real();
        cov_ok = arma::inv_sympd(proposed_cov_comb_inv.slice(b), proposed_cov_comb.slice(b));
        if(!cov_ok) {
          break;
        }
      }
      if(!cov_ok) {
        continue;
      }

      // The proposed model score is updated by the posterior kernel with the
      // proposed covariance, it's log determinant, inverse and the new 
      // likelihood covariance
      
      // Rcpp::Rcout << "\n\nSIGMA MH: Calculate posterior kernel.\n";
      
      proposed_model_score += sigmaLogKernel(k, 
                                             proposed_cov_comb_log_det,
                                             sigma_proposed,
                                             proposed_cov,
                                             proposed_cov_inv,
                                             proposed_cov_comb_inv
      );
      
      // The current model score is based on the posterior kernel for the current
      // covariance matrix
      
      current_model_score += sigmaLogKernel(k, 
                                        cov_comb_log_det.row(k).t(),
                                        sigma.col(k),
                                        cov.slice(k),
                                        cov_inv.slice(k),
                                        cov_comb_inv.slices(k * B + B_inds)
      );
      
      // Rcpp::Rcout << "\n\nSIGMA MH: MH step.\n";
      
      // Accept or reject
      u = arma::randu();
      acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));
      
      // Rcpp::Rcout << "\n\nSIGMA MH: acceptance probability calculated.";
      
    }
    if( (u < acceptance_prob) || (N_k(k) == 0) ){
      
      // Rcpp::Rcout << "SIGMA MH: Update count.\n";
      sigma_count(k)++;
      
      // Rcpp::Rcout << "SIGMA MH: Update sigma.\n";
      sigma.col(k) = sigma_proposed;
      Sigma_mat.slice(k).diag() = sigma_proposed;
      
      
      // Rcpp::Rcout << "SIGMA MH: Update cov.\n";
      cov.slice(k) = proposed_cov;
      cov_inv.slice(k) = proposed_cov_inv;
      cov_log_det(k) = proposed_cov_log_det;
      
      // Rcpp::Rcout << "SIGMA MH: Update combination items.\n";
      for(arma::uword b = 0; b < B; b++) {
        
        // Rcpp::Rcout << "SIGMA MH: Update combination covariance.\n";
        cov_comb.slice(k * B + b) = proposed_cov_comb.slice(b);
        // Rcpp::Rcout << "SIGMA MH: Update combination covariance determinant.\n";
        cov_comb_log_det(k, b) = proposed_cov_comb_log_det(b);
        
        // Rcpp::Rcout << "SIGMA MH: Update combination inverse.\n";
        cov_comb_inv.slice(k * B + b) = proposed_cov_comb_inv.slice(b);
      }
      // Rcpp::Rcout << "SIGMA MH: Updated combination items.\n";
    }
  }
};
void mvnSamplerSeparationStrategy::clusterMeanMetropolis() {
  
  double u = 0.0, proposed_model_score = 0.0, acceptance_prob = 0.0, current_model_score = 0.0;
  arma::vec mu_proposed(P);
  arma::mat proposed_mean_sum(P, B);
  
  mu_proposed.zeros();
  proposed_mean_sum.zeros();
  
  for(arma::uword k = 0; k < K ; k++) {
    if(N_k(k) == 0){
      mu_proposed = arma::mvnrnd(mu_0, (1.0/kappa) * cov.slice(k), 1);
      for(arma::uword b = 0; b < B; b++) {
        proposed_mean_sum.col(b) = mu_proposed + m.col(b);
        if(include_interaction) {
          proposed_mean_sum.col(b) += gamma.slice(b).col(k);
        }
      }
    } else {
      for(arma::uword p = 0; p < P; p++){
        // The proposal window is now a diagonal matrix of common entries.
        mu_proposed(p) = (arma::randn() * mu_proposal_window) + mu(p, k);
      }
      for(arma::uword b = 0; b < B; b++) {
        proposed_mean_sum.col(b) = mu_proposed + m.col(b);
        if(include_interaction) {
          proposed_mean_sum.col(b) += gamma.slice(b).col(k);
        }
      }

      // The prior is included in the kernel
      proposed_model_score = muLogKernel(k, mu_proposed, proposed_mean_sum);
      current_model_score = muLogKernel(k, mu.col(k), mean_sum.cols(k * B + B_inds));
      
      u = arma::randu();
      acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));
      
    }
    
    if((u < acceptance_prob) || (N_k(k) == 0)) {
      mu.col(k) = mu_proposed;
      mu_count(k)++;
      
      for(arma::uword b = 0; b < B; b++) {
        mean_sum.col(k * B + b) = proposed_mean_sum.col(b);
      }
      
    }
  }
};

// Block Metropolis-Hastings update for gamma(p, ., .) - see
// mvnSampler::interactionMetropolis() for the full derivation (identical
// here: same groupLikelihood() signature/semantics as mvnSampler).
void mvnSamplerSeparationStrategy::interactionMetropolis() {

  arma::uvec combined_group = labels * B + batch_vec;
  arma::uvec all_items = arma::regspace<arma::uvec>(0, N - 1);

  arma::vec cov_det_flat(K * B);
  for(arma::uword k = 0; k < K; k++) {
    for(arma::uword b = 0; b < B; b++) {
      cov_det_flat(k * B + b) = cov_comb_log_det(k, b);
    }
  }

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

    double current_model_score = groupLikelihood(all_items, combined_group, cov_det_flat, mean_sum, cov_comb_inv) + prior_current;
    double proposed_model_score = groupLikelihood(all_items, combined_group, cov_det_flat, proposed_mean_sum, cov_comb_inv) + prior_proposed;

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

void mvnSamplerSeparationStrategy::updateBatchCorrectedData() {

  arma::mat mu_mat = mu.cols(labels);
  arma::mat location_correction = m.cols(batch_vec);

  // See mvnSampler::updateBatchCorrectedData() for why gamma must be
  // subtracted here too whenever it's in use.
  if(include_interaction) {
    for(arma::uword n = 0; n < N; n++) {
      location_correction.col(n) += gamma.slice(batch_vec(n)).col(labels(n));
    }
  }

  Y = ((X_t - mu_mat - location_correction) / sqrt(S.cols(batch_vec)) + mu_mat).t();
}

void mvnSamplerSeparationStrategy::metropolisStep() {
  
  clusterMeanMetropolis();
  
  // Metropolis step for cluster parameters
  // Rcpp::Rcout << "\nR MH.";
  rMHStep();
  
  // Rcpp::Rcout << "\nSigma MH.";
  sigmaMHStep();
  
  // Rcpp::Rcout << "\nNew MH moves complete.";
 
  
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




//' @title Diagnostic R/sigma-only chain (internal)
//' @description Runs the R (correlation) and sigma (marginal SD)
//' Metropolis-Hastings steps of \code{mvnSamplerSeparationStrategy} in
//' isolation, for a single cluster/batch, to check their mixing/recovery
//' independently of the rest of the sampler. Not part of the public
//' workflow - see \code{\link{batchSemiSupervisedMixtureModel}} with
//' \code{type = "MVN_LKJ"} for the real model.
//' @param X The data matrix (items in rows).
//' @param rho_true The true correlation used only for reporting/comparison
//' by the caller; has no effect on the chain itself.
//' @param n_iter Number of iterations to run.
//' @param r_pw,sigma_pw Proposal windows for the R and sigma updates.
//' @return A named list with the sampled \code{r}/\code{sigma1} traces and
//' their acceptance rates.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List diagRSigmaOnlyChain2(arma::mat X, double rho_true, arma::uword n_iter, double r_pw, double sigma_pw) {
  uword K = 1, B = 1, P = 2;
  uword N = X.n_rows;
  uvec labels(N, fill::zeros);
  uvec batch_vec(N, fill::zeros);
  uvec fixed(N, fill::zeros);
  vec concentration = {1.0};

  mvnSamplerSeparationStrategy s(K, B, 0.3, r_pw, sigma_pw, 0.2, 40.0, labels, batch_vec, concentration, X, fixed, 0.01, 3.0, 1.0, false, 1.0);

  s.mu.col(0) = zeros<vec>(P);
  s.R.slice(0) = eye<mat>(P, P);
  s.R.slice(0)(0, 1) = rho_true;
  s.R.slice(0)(1, 0) = rho_true;
  s.r_log_det(0) = arma::log_det(s.R.slice(0)).real();
  s.sigma.col(0) = ones<vec>(P);
  s.Sigma_mat.slice(0) = eye<mat>(P, P);
  s.S.col(0) = ones<vec>(P) * 1.0;
  s.m.col(0) = zeros<vec>(P);
  s.cov.slice(0) = s.Sigma_mat.slice(0) * s.R.slice(0) * s.Sigma_mat.slice(0);
  s.matrixCombinations();
  s.updateWeights();

  vec r_trace(n_iter), sigma1_trace(n_iter);
  for(uword it = 0; it < n_iter; it++) {
    s.rMHStep();
    s.sigmaMHStep();
    r_trace(it) = s.R.slice(0)(0, 1);
    sigma1_trace(it) = s.sigma(0, 0);
  }
  return Rcpp::List::create(
    Rcpp::Named("r") = r_trace,
    Rcpp::Named("sigma1") = sigma1_trace,
    Rcpp::Named("r_acc") = (double) s.r_count(0) / n_iter,
    Rcpp::Named("sigma_acc") = (double) s.sigma_count(0) / n_iter
  );
}
