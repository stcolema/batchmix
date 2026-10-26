// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
#include "pdfs.cpp"
#include "genericFunctions.cpp"
#include "sampler.cpp"
#include "mvnSamplerSeparationStrategy.cpp"
#include "mvnSamplerMixed.cpp"
using namespace arma;

static const double LOD = -0.7;

static void setHyper(mvnSamplerSeparationStrategy& s, double kappa, double lam2, double eta, double xi_sd) {
  s.kappa = kappa; s.mu_0 = zeros<vec>(s.P); s.eta = eta; s.xi = xi_sd;
  s.delta_2 = 1.0; s.lambda_2 = lam2;
  s.batch_shift_prior_precision = 1.0 / (s.delta_2 * s.lambda_2);
}

static void simulate(mvnSamplerSeparationStrategy& s, mvnSamplerMixed* mx) {
  for (uword n = 0; n < s.N; n++) {
    uword k = s.labels(n), b = s.batch_vec(n);
    mat C = s.cov.slice(k);
    for (uword p = 0; p < s.P; p++) C(p, p) *= s.S(p, b);
    mat L = chol(C, "lower");
    vec x = s.mu.col(k) + s.m.col(b) + L * randn<vec>(s.P);
    if (mx == nullptr) { s.X.row(n) = x.t(); continue; }
    for (uword p = 0; p < s.P; p++) {
      if (mx->column_type(p) == 1) { mx->X_raw(n, p) = x(p) > 0 ? 1.0 : 0.0; }
      else if (x(p) < LOD) { mx->X_raw(n, p) = LOD; mx->censor_code(n, p) = 1; }
      else { mx->X_raw(n, p) = x(p); mx->censor_code(n, p) = 0; }
    }
    s.X.row(n) = x.t();   // working copy; overwritten by updateLatentData
  }
  if (mx != nullptr) mx->X_raw_t = mx->X_raw.t();
  s.X_t = s.X.t();
}

static vec summaries(mvnSamplerSeparationStrategy& s) {
  vec g(11);
  mat C0 = s.cov.slice(0);
  g(0) = s.mu(0, 0); g(1) = s.mu(1, 1);
  g(2) = std::log(det(C0)); g(3) = C0(0, 0);
  g(4) = s.R(0, 1, 0);
  g(5) = s.m(0, 0); g(6) = s.m(1, 1);
  g(7) = s.S(0, 0); g(8) = s.S(1, 1);
  g(9) = s.w(0);
  g(10) = s.sigma(0, 1);
  return g;
}

static void drawLabels(mvnSamplerSeparationStrategy& s) {
  vec w = randg<vec>(s.K, distr_param(1.0, 1.0)); w /= accu(w); s.w = w;
  for (uword n = 0; n < s.N; n++) s.labels(n) = (randu() < w(0)) ? 0 : 1;
}

// [[Rcpp::export]]
Rcpp::List gewekeSep(bool mixed, unsigned n_mc, unsigned n_sc, unsigned thin,
                     double r_pw, double sigma_pw, double mu_pw, double m_pw, double S_pw, double eta) {
  uword K = 2, B = 2, P = 2, N = 10;
  uvec batch_vec(N), fixed(N, fill::zeros), labels(N, fill::zeros);
  for (uword n = 0; n < N; n++) batch_vec(n) = n % B;
  vec conc = {1.0, 1.0};
  mat X = randn<mat>(N, P);
  uvec ctype = {0, 1};
  umat ccode(N, P, fill::zeros);
  mvnSamplerSeparationStrategy* sp; mvnSamplerMixed* mx = nullptr;
  if (mixed) {
    mx = new mvnSamplerMixed(K, B, mu_pw, r_pw, sigma_pw, m_pw, S_pw, labels, batch_vec, conc, X, fixed, 0.5, 3.0, 1.0, false, eta, ctype, ccode, false);
    sp = mx;
  } else {
    sp = new mvnSamplerSeparationStrategy(K, B, mu_pw, r_pw, sigma_pw, m_pw, S_pw, labels, batch_vec, conc, X, fixed, 0.5, 3.0, 1.0, false, eta, false);
  }
  mvnSamplerSeparationStrategy& s = *sp;
  setHyper(s, 1.0, 0.5, eta, 0.7);
  mat mc(n_mc, 11), sc(n_sc, 11);
  for (unsigned i = 0; i < n_mc; i++) {
    s.sampleFromPriors(); s.matrixCombinations(); drawLabels(s);
    mc.row(i) = summaries(s).t();
  }
  s.sampleFromPriors(); s.matrixCombinations(); drawLabels(s);
  unsigned total = 0;
  for (unsigned i = 0; i < n_sc; i++) {
    for (unsigned t = 0; t < thin; t++) {
      simulate(s, mx);
      s.updateLatentData();
      s.updateWeights();
      s.metropolisStep();
      s.updateAllocation();
      total++;
    }
    sc.row(i) = summaries(s).t();
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("mc") = mc, Rcpp::Named("sc") = sc,
    Rcpp::Named("r_acc") = conv_to<vec>::from(s.r_count) / total,
    Rcpp::Named("sigma_acc") = conv_to<vec>::from(s.sigma_count) / total,
    Rcpp::Named("mu_acc") = conv_to<vec>::from(s.mu_count) / total,
    Rcpp::Named("m_acc") = conv_to<vec>::from(s.m_count) / total,
    Rcpp::Named("S_acc") = conv_to<vec>::from(s.S_count) / total);
  delete sp;
  return out;
}
