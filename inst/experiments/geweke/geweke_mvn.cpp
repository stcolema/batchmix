// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
#include "pdfs.cpp"
#include "genericFunctions.cpp"
#include "sampler.cpp"
#include "mvnSampler.cpp"
#include "mvtSampler.cpp"
using namespace arma;

static void setHyper(mvnSampler& s, double kappa, double lam2) {
  s.kappa = kappa;
  s.xi = zeros<vec>(s.P);
  s.scale = eye<mat>(s.P, s.P) * 0.5;
  s.nu = s.P + 2;
  s.delta_2 = 1.0;
  s.lambda_2 = lam2;
  s.batch_shift_prior_precision = 1.0 / (s.delta_2 * s.lambda_2);
}

static void simulateData(mvnSampler& s, const vec& tdf, bool is_t) {
  // independent reconstruction from raw parameters (not the sampler's caches)
  for (uword n = 0; n < s.N; n++) {
    uword k = s.labels(n), b = s.batch_vec(n);
    mat C = s.cov.slice(k);
    for (uword p = 0; p < s.P; p++) C(p, p) *= s.S(p, b);
    mat L = chol(C, "lower");
    vec z = randn<vec>(s.P);
    double scale = 1.0;
    if (is_t) scale = 1.0 / std::sqrt(randg<double>(distr_param(0.5 * tdf(k), 2.0 / tdf(k))));
    vec x = s.mu.col(k) + s.m.col(b) + scale * (L * z); if (s.include_interaction) x += s.gamma.slice(b).col(k);
    s.X.row(n) = x.t();
  }
  s.X_t = s.X.t();
}

static vec summaries(mvnSampler& s, const vec& tdf, bool is_t) {
  mat C0 = s.cov.slice(0);
  vec g(12 + (is_t ? 1 : 0));
  g(0) = s.mu(0, 0); g(1) = s.mu(1, 1);
  g(2) = std::log(det(C0)); g(3) = C0(0, 0);
  g(4) = C0(0, 1) / std::sqrt(C0(0, 0) * C0(1, 1));
  g(5) = s.m(0, 0); g(6) = s.m(1, 1);
  g(7) = s.S(0, 0); g(8) = s.S(1, 1); 
  g(9) = (double) accu(s.labels == 0);
  g(10) = s.w(0);
  g(11) = s.include_interaction ? s.gamma(0, 0, 0) : s.mu(0, 0) * s.m(0, 0);
  if (s.include_interaction) g(10) = std::log(s.tau2_interaction(0));
  if (is_t) g(12) = tdf(0);
  return g;
}

// [[Rcpp::export]]
Rcpp::List gewekeMVN(bool is_t, unsigned n_mc, unsigned n_sc, unsigned thin,
                     double mu_pw, double cov_pw, double m_pw, double S_pw, double df_pw,
                     bool sample_m_scale, bool sample_s_scale, unsigned seed, bool inter, bool miss) {
  arma::arma_rng::set_seed(seed);
  uword K = 2, B = 2, P = 2, N = 10;
  uvec batch_vec(N), fixed(N, fill::zeros), labels(N, fill::zeros);
  for (uword n = 0; n < N; n++) batch_vec(n) = n % B;
  vec conc = {1.0, 1.0};
  mat X = randn<mat>(N, P);
  mvtSampler* ts = nullptr;
  mvnSampler* sp = nullptr;
  if (is_t) { ts = new mvtSampler(K, B, mu_pw, cov_pw, m_pw, S_pw, df_pw, labels, batch_vec, conc, X, fixed, 0.5, 3.0, 1.0, sample_m_scale, sample_s_scale); sp = ts; }
  else sp = new mvnSampler(K, B, mu_pw, cov_pw, m_pw, S_pw, labels, batch_vec, conc, X, fixed, 0.5, 3.0, 1.0, sample_m_scale, sample_s_scale);
  mvnSampler& s = *sp;
  setHyper(s, 1.0, 0.5);
  s.s_scale_proposal_window = 0.5;
  s.initialiseInteraction(inter, 0.3, 3.0, 2.0);
  if (miss) { std::vector<uword> it; for (uword n = 0; n < N; n++) if (n % 3 == 0) it.push_back(n); s.items_to_augment = conv_to<uvec>::from(it); for (uword n : it) s.X_raw_t(1, n) = datum::nan; }
  vec dummy_df = ones<vec>(K) * 10.0;
  auto tdf = [&]() -> vec { return is_t ? ts->t_df : dummy_df; };

  unsigned nsum = 12 + (is_t ? 1 : 0);
  mat mc(n_mc, nsum), sc(n_sc, nsum);

  // marginal-conditional: independent prior draws
  for (unsigned i = 0; i < n_mc; i++) {
    if (sample_s_scale) { s.rho = 2.0 + randg<double>(distr_param(s.a_s, 1.0 / s.b_s)); s.theta = s.s_scale_prior_mean * (s.rho - 1.0); }
    if (inter) for (uword p = 0; p < P; p++) s.tau2_interaction(p) = rInvGamma(s.a_gamma, s.b_gamma);
    s.sampleFromPriors();
    s.matrixCombinations();
    vec w = randg<vec>(K, distr_param(1.0, 1.0)); w /= accu(w);
    s.w = w;
    for (uword n = 0; n < N; n++) s.labels(n) = (randu() < w(0)) ? 0 : 1;
    mc.row(i) = summaries(s, tdf(), is_t).t();
  }

  // successive-conditional
  if (sample_s_scale) { s.rho = 2.0 + randg<double>(distr_param(s.a_s, 1.0 / s.b_s)); s.theta = s.s_scale_prior_mean * (s.rho - 1.0); }
  if (inter) for (uword p = 0; p < P; p++) s.tau2_interaction(p) = rInvGamma(s.a_gamma, s.b_gamma);
  s.sampleFromPriors();
  s.matrixCombinations();
  { vec w = randg<vec>(K, distr_param(1.0, 1.0)); w /= accu(w); s.w = w;
    for (uword n = 0; n < N; n++) s.labels(n) = (randu() < w(0)) ? 0 : 1; }
  double acc_cov = 0, acc_mu = 0, acc_m = 0, acc_S = 0;
  unsigned total = 0;
  for (unsigned i = 0; i < n_sc; i++) {
    for (unsigned t = 0; t < thin; t++) {
      simulateData(s, tdf(), is_t);
      s.updateLatentData();
      s.updateWeights();
      s.metropolisStep();
      s.updateAllocation();
      total++;
    }
    sc.row(i) = summaries(s, tdf(), is_t).t();
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("mc") = mc, Rcpp::Named("sc") = sc,
    Rcpp::Named("cov_acc") = conv_to<vec>::from(s.cov_count) / total,
    Rcpp::Named("mu_acc") = conv_to<vec>::from(s.mu_count) / total,
    Rcpp::Named("m_acc") = conv_to<vec>::from(s.m_count) / total,
    Rcpp::Named("S_acc") = conv_to<vec>::from(s.S_count) / total);
  if (is_t) out["df_acc"] = conv_to<vec>::from(ts->t_df_count) / total;
  delete sp;
  return out;
}
