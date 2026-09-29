// priorOnlyChains.cpp
// =============================================================================
// Test-only diagnostic chains. Each subclasses a concrete sampler and replaces
// the data log-likelihood (groupLikelihood()) with zero, so every Metropolis
// step targets the PRIOR alone. Comparing the resulting marginals with direct
// prior draws ("prior recovery"; e.g. Talts et al., 2018, arXiv:1804.06788,
// Sec. 2 on checking a sampler against its own prior) exposes errors in a
// kernel's prior terms, Jacobians and empty-cluster handling that a
// data-driven test cannot separate from the likelihood.
# include <RcppArmadillo.h>
# include "sampler.h"
# include "mvnSamplerSeparationStrategy.h"

using namespace Rcpp ;
using namespace arma ;

namespace {

class priorOnlyLKJSampler : public mvnSamplerSeparationStrategy {
public:
  using mvnSamplerSeparationStrategy::mvnSamplerSeparationStrategy;

  double groupLikelihood(arma::uvec, arma::uvec, arma::vec, arma::mat, arma::cube) override {
    return 0.0;
  }
};

}

//' @title Prior-only MVN_LKJ chain (internal)
//' @description Runs the MVN_LKJ Metropolis steps with the data likelihood
//' switched off, so the chain's stationary distribution is the prior. Items
//' in \code{X} only serve to make clusters occupied (those in
//' \code{labels}) or empty (every other cluster).
//' @param X Data matrix (only its dimensions and empirical-Bayes prior
//' centre matter).
//' @param K,B Number of clusters and batches.
//' @param labels Cluster label (0-indexed) of each row of \code{X}.
//' @param batch_vec Batch label (0-indexed) of each row of \code{X}.
//' @param n_iter Number of sweeps.
//' @param eta LKJ concentration.
//' @param r_pw,sigma_pw,mu_pw,m_pw,S_pw Proposal windows, on the scale the
//' sampler uses internally.
//' @return A list of traces, acceptance counts and initial-state
//' diagnostics.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List priorOnlyLKJChain(
  arma::mat X,
  arma::uword K,
  arma::uword B,
  arma::uvec labels,
  arma::uvec batch_vec,
  arma::uword n_iter,
  double eta,
  double r_pw,
  double sigma_pw,
  double mu_pw,
  double m_pw,
  double S_pw
) {
  uword N = X.n_rows, P = X.n_cols;
  uvec fixed(N, fill::zeros);
  vec concentration = ones<vec>(K);

  priorOnlyLKJSampler s(K, B, mu_pw, r_pw, sigma_pw, m_pw, S_pw, labels,
    batch_vec, concentration, X, fixed, 0.01, 3.0, 1.0, false, eta);

  s.sampleFromPriors();
  s.matrixCombinations();

  // State straight after initialisation, before any Metropolis step.
  vec r_log_det_init = s.r_log_det;
  vec r_log_det_true(K);
  for(uword k = 0; k < K; k++) {
    r_log_det_true(k) = arma::log_det(s.R.slice(k)).real();
  }

  cube sigma_trace(P, K, n_iter), mu_trace(P, K, n_iter),
    S_trace(P, B, n_iter), cov_trace(P, P * K, n_iter);
  mat r_trace(n_iter, K);

  for(uword it = 0; it < n_iter; it++) {
    s.updateWeights();
    s.metropolisStep();

    sigma_trace.slice(it) = s.sigma;
    mu_trace.slice(it) = s.mu;
    S_trace.slice(it) = s.S;
    for(uword k = 0; k < K; k++) {
      cov_trace.slice(it).cols(k * P, (k + 1) * P - 1) = s.cov.slice(k);
      r_trace(it, k) = s.R.slice(k)(0, 1);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("sigma") = sigma_trace,
    Rcpp::Named("mu") = mu_trace,
    Rcpp::Named("S") = S_trace,
    Rcpp::Named("cov") = cov_trace,
    Rcpp::Named("r") = r_trace,
    Rcpp::Named("mu_0") = s.mu_0,
    Rcpp::Named("kappa") = s.kappa,
    Rcpp::Named("beta") = s.beta,
    Rcpp::Named("xi") = s.xi,
    Rcpp::Named("r_log_det_init") = r_log_det_init,
    Rcpp::Named("r_log_det_true") = r_log_det_true,
    Rcpp::Named("mu_count") = conv_to<vec>::from(s.mu_count),
    Rcpp::Named("r_count") = conv_to<vec>::from(s.r_count),
    Rcpp::Named("sigma_count") = conv_to<vec>::from(s.sigma_count)
  );
}

//' @title Partial-pooling batch-weight chain with fixed hyperparameters (internal)
//' @description Runs \code{updatePartialPoolingWeights()} repeatedly,
//' resetting the population mean and variance to fixed values before every
//' sweep so that the ALR entries of each batch are updated by the Metropolis
//' kernel alone, targeting their exact conditional
//' \eqn{p(\eta_b \mid counts_b) \propto \exp(c_b \eta_b - N_b \log(1 + e^{\eta_b}))
//' N(\eta_b; \mu, \tau^2)} (for K = 2). A batch with no items follows its
//' prior N(mu, tau2).
//' @param labels,batch_vec 0-indexed cluster and batch labels of the items.
//' @param K,B Number of clusters and batches.
//' @param n_iter Number of sweeps.
//' @param mu,tau2 The fixed hyperparameters (for every free coordinate).
//' @param eta_pw The proposal window (sd) of each scalar random walk.
//' @return A list with the \code{eta} trace (B x (K - 1) x n_iter), the raw
//' \code{eta_count}, and \code{eta_moves_per_sweep}.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List partialPoolingWeightChain(
  arma::uvec labels,
  arma::uvec batch_vec,
  arma::uword K,
  arma::uword B,
  arma::uword n_iter,
  double mu,
  double tau2,
  double eta_pw
) {
  uword N = labels.n_elem;
  arma::mat X = arma::randn<arma::mat>(N, 2);
  uvec fixed(N, fill::zeros);
  vec concentration = ones<vec>(K);

  priorOnlyLKJSampler s(K, B, 0.3, 1.0, 5.0, 0.3, 5.0, labels, batch_vec,
    concentration, X, fixed, 0.01, 3.0, 1.0, false, 1.0);
  s.initialiseBatchWeightPrior(1, arma::regspace<arma::vec>(0, B - 1), 1.0, 1.0,
    eta_pw, false, 0.1, 2.0, 1.0, 10.0);

  cube eta_trace(B, K - 1, n_iter);
  for(uword it = 0; it < n_iter; it++) {
    s.pp_mu.fill(mu);
    s.pp_tau2.fill(tau2);
    s.updateWeights();
    eta_trace.slice(it) = s.eta_alr;
  }

  return Rcpp::List::create(
    Rcpp::Named("eta") = eta_trace,
    Rcpp::Named("eta_count") = conv_to<vec>::from(s.eta_count),
    Rcpp::Named("eta_moves_per_sweep") = s.eta_moves_per_sweep
  );
}
