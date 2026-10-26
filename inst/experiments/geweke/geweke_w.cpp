// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
#include "pdfs.cpp"
#include "genericFunctions.cpp"
#include "sampler.cpp"
#include "mvnSampler.cpp"
using namespace arma;

// draw the weight-block state from its prior
static void drawPrior(mvnSampler& s, uword type, uword kernel, bool hyper) {
  uword B = s.B, K = s.K;
  if (type == 1) {
    vec mu(K); for (uword k = 0; k < K; k++) mu(k) = randn() * s.pp_mu_prior_sd;
    mu -= mean(mu); s.pp_mu = mu;
    for (uword k = 0; k < K; k++) s.pp_tau2(k) = rInvGamma(s.pp_tau2_prior_shape, s.pp_tau2_prior_rate);
    for (uword b = 0; b < B; b++) for (uword k = 0; k < K; k++) s.eta_logit(b, k) = s.pp_mu(k) + std::sqrt(s.pp_tau2(k)) * randn();
  } else {
    if (hyper) {
      s.gp_tau2 = rInvGamma(s.gp_tau2_prior_shape, s.gp_tau2_prior_rate);
      if (kernel == 0) s.gp_length_scale = rInvGamma(s.gp_length_scale_prior_shape, s.gp_length_scale_prior_rate);
    }
    s.buildWellConditionedGPChol(s.gp_tau2, s.gp_length_scale, s.gp_cov, s.gp_chol, s.gp_jitter);
    if (kernel == 0) {
      vec beta(K); for (uword k = 0; k < K; k++) beta(k) = randn() * s.pp_mu_prior_sd;
      beta -= mean(beta); s.gp_beta = beta;
    } else s.gp_beta.zeros();
    mat L = trimatl(s.gp_chol);
    for (uword k = 0; k < K; k++) s.eta_logit.col(k) = s.gp_beta(k) + L * randn<vec>(B);
  }
  s.updateSimplexFromLogits();
}

static void drawLabels(mvnSampler& s) {
  for (uword n = 0; n < s.N; n++) {
    uword b = s.batch_vec(n);
    double u = randu(), c = 0; uword k = 0;
    for (k = 0; k < s.K - 1; k++) { c += s.w_batch(b, k); if (u < c) break; }
    s.labels(n) = k;
  }
}

static vec summ(mvnSampler& s, uword type, uword kernel) {
  vec g(8);
  g(0) = s.eta_logit(0, 0) - s.eta_logit(0, 1);
  g(1) = s.eta_logit(1, 0) - s.eta_logit(1, 2);
  g(2) = s.eta_logit(s.B - 1, 1) - s.eta_logit(s.B - 1, 0);
  g(3) = (type == 1) ? s.pp_mu(0) : s.gp_beta(0);
  g(4) = (type == 1) ? std::log(s.pp_tau2(0)) : std::log(s.gp_tau2);
  g(5) = (type == 1) ? std::log(s.pp_tau2(1)) : std::log(s.gp_length_scale);
  g(6) = s.w_batch(0, 0);
  g(7) = s.w_batch(s.B - 1, 2);
  return g;
}

// [[Rcpp::export]]
Rcpp::List gewekeW(unsigned type, unsigned kernel, bool hyper, unsigned n_mc, unsigned n_sc, unsigned thin,
                   double eta_pw, double hyper_pw, double pp_mu_sd, unsigned per_batch) {
  uword K = 3, B = 5, N = B * per_batch;
  uvec batch_vec(N), fixed(N, fill::zeros), labels(N, fill::zeros);
  for (uword n = 0; n < N; n++) batch_vec(n) = n % B;
  vec conc = ones<vec>(K);
  mat X = randn<mat>(N, 2);
  mvnSampler s(K, B, 0.5, 30, 0.5, 30, labels, batch_vec, conc, X, fixed, 0.5, 3.0, 1.0, false, false);
  vec coords = {0.0, 0.15, 0.4, 0.7, 1.0};
  s.initialiseBatchWeightPrior(type, coords, 1.0, 0.5, eta_pw, hyper, hyper_pw, 2.0, 1.0, pp_mu_sd, kernel);
  // initialiseBatchWeightPrior sets N_b from the constructor (labels irrelevant)
  mat mc(n_mc, 8), sc(n_sc, 8);
  for (unsigned i = 0; i < n_mc; i++) { drawPrior(s, type, kernel, hyper); mc.row(i) = summ(s, type, kernel).t(); }
  drawPrior(s, type, kernel, hyper);
  double eta_acc_den = 0; 
  for (unsigned i = 0; i < n_sc; i++) {
    for (unsigned t = 0; t < thin; t++) { drawLabels(s); s.updateWeights(); eta_acc_den++; }
    sc.row(i) = summ(s, type, kernel).t();
  }
  return Rcpp::List::create(Rcpp::Named("mc") = mc, Rcpp::Named("sc") = sc,
    Rcpp::Named("eta_acc") = conv_to<vec>::from(s.eta_count) / (eta_acc_den * n_sc * 0 + (double) n_sc * thin * s.eta_moves_per_sweep),
    Rcpp::Named("hyper_acc") = (double) s.gp_hyperparameter_count / ((double) n_sc * thin));
}
