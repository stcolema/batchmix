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
# include "mvnSampler.h"
# include "mvtSampler.h"
# include <memory>

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

// The plain MVN and MVT samplers with the data likelihood replaced by zero.
// Every data term of every Metropolis kernel (and of the label update) goes
// through the members overridden here, so the kernels' prior terms, Jacobians
// and proposal corrections are what is left to be tested.
class priorOnlyMVNSampler : public mvnSampler {
public:
  using mvnSampler::mvnSampler;

  double groupLikelihood(arma::uvec, arma::uvec, arma::vec, arma::mat, arma::cube) override {
    return 0.0;
  }

  arma::vec itemLogLikelihood(arma::vec, arma::uword) override {
    return arma::zeros<arma::vec>(K);
  }
};

class priorOnlyMVTSampler : public mvtSampler {
public:
  using mvtSampler::mvtSampler;

  double clusterLikelihood(double, arma::uvec, arma::vec, arma::mat, arma::cube) override {
    return 0.0;
  }

  double batchLikelihood(arma::uvec, arma::uvec, arma::vec, arma::vec, arma::mat, arma::cube) override {
    return 0.0;
  }

  double dfLogLikelihood(arma::uword, double, double) override {
    return 0.0;
  }

  arma::vec itemLogLikelihood(arma::vec, arma::uword) override {
    return arma::zeros<arma::vec>(K);
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
//' resetting the population mean and variances to fixed values before every
//' sweep so that each batch's logits are updated by the Metropolis and shift
//' Gibbs kernels alone, targeting their exact conditional. The weights depend
//' on the logits only through differences, so for two classes the logit
//' difference \eqn{d = \eta_1 - \eta_2} of a batch with counts \eqn{c_1} and
//' \eqn{N} has density
//' \eqn{\propto \exp(c_1 d - N \log(1 + e^d)) N(d; \mu_1 - \mu_2, \tau_1^2 +
//' \tau_2^2)}; a batch with no items follows the prior.
//' @param labels,batch_vec 0-indexed cluster and batch labels of the items.
//' @param K,B Number of clusters and batches.
//' @param n_iter Number of sweeps.
//' @param mu,tau2 K-vectors: the fixed population mean and class variances.
//' @param eta_pw The proposal window (sd) of each scalar random walk.
//' @return A list with the \code{eta} trace (B x K x n_iter), the raw
//' \code{eta_count}, and \code{eta_moves_per_sweep}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List partialPoolingWeightChain(
  arma::uvec labels,
  arma::uvec batch_vec,
  arma::uword K,
  arma::uword B,
  arma::uword n_iter,
  arma::vec mu,
  arma::vec tau2,
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

  cube eta_trace(B, K, n_iter);
  for(uword it = 0; it < n_iter; it++) {
    s.pp_mu = mu;
    s.pp_tau2 = tau2;
    s.updateWeights();
    eta_trace.slice(it) = s.eta_logit;
  }

  return Rcpp::List::create(
    Rcpp::Named("eta") = eta_trace,
    Rcpp::Named("eta_count") = conv_to<vec>::from(s.eta_count),
    Rcpp::Named("eta_moves_per_sweep") = s.eta_moves_per_sweep
  );
}

//' @title Test helper: GP batch-weight update with no data
//' @description Runs \code{updateWeights()} for the GP batch-weight prior with
//' no items, so the chain targets the prior: each class's logits follow
//' \eqn{N(0, \Sigma)} with \eqn{\Sigma} from \code{gpKernelMatrix()}
//' (the common per-batch shift is redrawn exactly). Used to check that the
//' kernels, including the non-stationary ones, are sampled correctly.
//' @param K,B Number of clusters and batches.
//' @param n_iter Number of sweeps.
//' @param coords B-vector of (already rescaled) batch coordinates.
//' @param kernel_type 0, 1 or 2; see \code{gpKernelMatrix()}.
//' @param tau2,length_scale Kernel parameters.
//' @param level_sd The prior sd of the diffuse level (and slope).
//' @param eta_pw Proposal window of the block random walk.
//' @param sample_hyperparameters Update \code{tau2} and \code{length_scale}
//' by Metropolis-Hastings, so that they follow their hyperpriors.
//' @param hyper_pw Proposal window of that update.
//' @return A list with the \code{eta} trace (B x K x n_iter),
//' \code{eta_count}, the traces of the hyperparameters and the shape and rate
//' of their priors.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List gpWeightChain(
  arma::uword K,
  arma::uword B,
  arma::uword n_iter,
  arma::vec coords,
  arma::uword kernel_type,
  double tau2,
  double length_scale,
  double level_sd,
  double eta_pw,
  bool sample_hyperparameters = false,
  double hyper_pw = 0.1
) {
  // The sampler's constructor needs at least one item; two dummy items are
  // dropped again below so that no batch carries any data.
  arma::uvec labels = arma::zeros<arma::uvec>(2), batch_vec = arma::zeros<arma::uvec>(2),
    fixed = arma::zeros<arma::uvec>(2);
  arma::mat X = arma::randn<arma::mat>(2, 2);
  vec concentration = ones<vec>(K);

  priorOnlyLKJSampler s(K, B, 0.3, 1.0, 5.0, 0.3, 5.0, labels, batch_vec,
    concentration, X, fixed, 0.01, 3.0, 1.0, false, 1.0);
  s.initialiseBatchWeightPrior(2, coords, tau2, length_scale,
    eta_pw, sample_hyperparameters, hyper_pw, 2.0, 1.0, level_sd, kernel_type);
  s.labels = arma::zeros<arma::uvec>(0);
  s.batch_vec = arma::zeros<arma::uvec>(0);
  s.N_b = arma::zeros<arma::uvec>(B);
  s.members = arma::zeros<arma::umat>(0, K);

  cube eta_trace(B, K, n_iter);
  vec tau2_trace(n_iter), length_trace(n_iter);
  for(uword it = 0; it < n_iter; it++) {
    s.updateWeights();
    eta_trace.slice(it) = s.eta_logit;
    tau2_trace(it) = s.gp_tau2;
    length_trace(it) = s.gp_length_scale;
  }

  return Rcpp::List::create(
    Rcpp::Named("eta") = eta_trace,
    Rcpp::Named("eta_count") = conv_to<vec>::from(s.eta_count),
    Rcpp::Named("gp_tau2") = tau2_trace,
    Rcpp::Named("gp_length_scale") = length_trace,
    Rcpp::Named("gp_tau2_prior_shape") = s.gp_tau2_prior_shape,
    Rcpp::Named("gp_tau2_prior_rate") = s.gp_tau2_prior_rate,
    Rcpp::Named("gp_length_scale_prior_shape") = s.gp_length_scale_prior_shape,
    Rcpp::Named("gp_length_scale_prior_rate") = s.gp_length_scale_prior_rate,
    Rcpp::Named("gp_hyperparameter_count") = s.gp_hyperparameter_count
  );
}

//' @title Prior-only chain for the plain MVN and MVT samplers (internal)
//' @description Runs the sampler's sweep (weights, Metropolis steps and, if
//' asked, the labels) with the data likelihood switched off, so the chain's
//' stationary distribution is the prior. The initial state is a draw from the
//' prior, so the cluster, batch and hyperparameter blocks are stationary
//' from the first sweep when the labels are held fixed.
//' @param type \code{"MVN"} or \code{"MVT"}.
//' @param X Data matrix (only its dimensions and the empirical-Bayes prior
//' centre and scale matter).
//' @param K,B Number of clusters and batches.
//' @param labels Cluster label (0-indexed) of each row of \code{X}.
//' @param batch_vec Batch label (0-indexed) of each row of \code{X}.
//' @param n_iter Number of sweeps.
//' @param mu_pw,cov_pw,m_pw,S_pw,t_df_pw,s_scale_pw Proposal windows, on the
//' scale the sampler uses internally.
//' @param m_scale,rho,theta The batch shift and scale hyperparameters.
//' @param sample_m_scale,sample_s_scale Sample the shift prior's variance and
//' the scale prior's concentration.
//' @param concentration K-vector, the Dirichlet concentration of the weights.
//' @param update_allocation Update the labels (they follow the prior) rather
//' than holding them fixed.
//' @return A list of traces and the prior's hyperparameters.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List priorOnlyChain(
  std::string type,
  arma::mat X,
  arma::uword K,
  arma::uword B,
  arma::uvec labels,
  arma::uvec batch_vec,
  arma::uword n_iter,
  double mu_pw,
  double cov_pw,
  double m_pw,
  double S_pw,
  double t_df_pw,
  double m_scale,
  double rho,
  double theta,
  bool sample_m_scale,
  bool sample_s_scale,
  double s_scale_pw,
  arma::vec concentration,
  bool update_allocation
) {
  const bool is_t = (type == "MVT");
  if(!is_t && type != "MVN") {
    Rcpp::stop("type must be 'MVN' or 'MVT'.");
  }
  uword N = X.n_rows, P = X.n_cols;
  uvec fixed(N, fill::zeros);

  std::unique_ptr<mvnSampler> sampler_ptr;
  if(is_t) {
    sampler_ptr.reset(new priorOnlyMVTSampler(K, B, mu_pw, cov_pw, m_pw, S_pw, t_df_pw,
      labels, batch_vec, concentration, X, fixed, m_scale, rho, theta,
      sample_m_scale, sample_s_scale));
  } else {
    sampler_ptr.reset(new priorOnlyMVNSampler(K, B, mu_pw, cov_pw, m_pw, S_pw,
      labels, batch_vec, concentration, X, fixed, m_scale, rho, theta,
      sample_m_scale, sample_s_scale));
  }
  mvnSampler& s = *sampler_ptr;
  mvtSampler* t_sampler = is_t ? dynamic_cast<mvtSampler*>(sampler_ptr.get()) : nullptr;
  s.s_scale_proposal_window = s_scale_pw;

  s.sampleFromPriors();
  s.matrixCombinations();

  cube mu_trace(P, K, n_iter), cov_trace(P, P * K, n_iter),
    m_trace(P, B, n_iter), S_trace(P, B, n_iter);
  vec lambda_2_trace(n_iter), rho_trace(n_iter);
  umat label_trace(n_iter, N);
  mat w_trace(n_iter, K), t_df_trace(n_iter, K, fill::zeros);

  for(uword it = 0; it < n_iter; it++) {
    s.updateWeights();
    s.metropolisStep();
    if(update_allocation) {
      s.updateAllocation();
    }

    mu_trace.slice(it) = s.mu;
    m_trace.slice(it) = s.m;
    S_trace.slice(it) = s.S;
    for(uword k = 0; k < K; k++) {
      cov_trace.slice(it).cols(k * P, (k + 1) * P - 1) = s.cov.slice(k);
    }
    lambda_2_trace(it) = s.lambda_2;
    rho_trace(it) = s.rho;
    label_trace.row(it) = s.labels.t();
    w_trace.row(it) = s.w.t();
    if(t_sampler != nullptr) {
      t_df_trace.row(it) = t_sampler->t_df.t();
    }
  }

  // List::create() takes at most 20 arguments in older Rcpp, so assign one by one
  Rcpp::List out;
  out["mu"] = mu_trace;
  out["cov"] = cov_trace;
  out["m"] = m_trace;
  out["S"] = S_trace;
  out["lambda_2"] = lambda_2_trace;
  out["rho"] = rho_trace;
  out["labels"] = label_trace;
  out["w"] = w_trace;
  out["t_df"] = t_df_trace;
  out["xi"] = s.xi;
  out["kappa"] = s.kappa;
  out["nu"] = s.nu;
  out["scale"] = s.scale;
  out["delta_2"] = s.delta_2;
  out["a"] = s.a;
  out["b"] = s.b;
  out["a_s"] = s.a_s;
  out["b_s"] = s.b_s;
  out["s_scale_prior_mean"] = s.s_scale_prior_mean;
  out["S_loc"] = s.S_loc;
  out["mu_count"] = conv_to<vec>::from(s.mu_count);
  out["cov_count"] = conv_to<vec>::from(s.cov_count);
  out["m_count"] = conv_to<vec>::from(s.m_count);
  out["S_count"] = conv_to<vec>::from(s.S_count);
  out["s_scale_count"] = (double) s.s_scale_count;
  if(t_sampler != nullptr) {
    out["psi"] = t_sampler->psi;
    out["chi"] = t_sampler->chi;
    out["t_loc"] = t_sampler->t_loc;
    out["t_df_count"] = conv_to<vec>::from(t_sampler->t_df_count);
  }
  return out;
}

//' @title Prior-only chain for the batch-by-cluster interaction (internal)
//' @description Runs the interaction term's Metropolis step and the Gibbs
//' update of its variance with no items, so that no data likelihood enters and
//' the chain targets the prior: each \code{gamma[p, , ]} is Gaussian on the
//' sum-to-zero subspace with variance \code{tau2_interaction[p]}, itself
//' \code{InvGamma(a_gamma, b_gamma)}.
//' @param type \code{"MVN"} or \code{"MVT"} (which has its own interaction step).
//' @param P,K,B Number of features, clusters and batches.
//' @param n_iter Number of sweeps.
//' @param gamma_pw Proposal window (sd) of the block random walk.
//' @param a_gamma,b_gamma Shape and rate of the Inverse-Gamma hyperprior.
//' @return A list with the traces of \code{tau2} (n_iter x P), of the sum of
//' squares of \code{gamma[p, , ]} (n_iter x P) and of \code{gamma[1, , ]}
//' (K x B x n_iter), and the acceptance counts.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List priorOnlyInteractionChain(
  std::string type,
  arma::uword P,
  arma::uword K,
  arma::uword B,
  arma::uword n_iter,
  double gamma_pw,
  double a_gamma,
  double b_gamma
) {
  const bool is_t = (type == "MVT");
  if(!is_t && type != "MVN") {
    Rcpp::stop("type must be 'MVN' or 'MVT'.");
  }
  // The sampler needs items to be constructed; they are removed below
  arma::mat X = arma::randn<arma::mat>(4, P);
  arma::uvec labels = arma::zeros<arma::uvec>(4), batch_vec = arma::zeros<arma::uvec>(4),
    fixed = arma::zeros<arma::uvec>(4);
  vec concentration = ones<vec>(K);

  std::unique_ptr<mvnSampler> sampler_ptr;
  if(is_t) {
    sampler_ptr.reset(new priorOnlyMVTSampler(K, B, 1.0, 10.0, 0.1, 5.0, 10.0,
      labels, batch_vec, concentration, X, fixed, 0.1, 3.0, 1.0, false, false));
  } else {
    sampler_ptr.reset(new priorOnlyMVNSampler(K, B, 1.0, 10.0, 0.1, 5.0,
      labels, batch_vec, concentration, X, fixed, 0.1, 3.0, 1.0, false, false));
  }
  mvnSampler& s = *sampler_ptr;
  s.initialiseInteraction(true, gamma_pw, a_gamma, b_gamma);
  s.sampleFromPriors();
  s.matrixCombinations();

  // No items: the data likelihood of every proposal is exactly zero
  s.labels = arma::zeros<arma::uvec>(0);
  s.batch_vec = arma::zeros<arma::uvec>(0);
  s.X_t = arma::zeros<arma::mat>(P, 0);
  s.N = 0;

  mat tau2_trace(n_iter, P), ss_trace(n_iter, P);
  cube gamma_trace(K, B, n_iter);
  for(uword it = 0; it < n_iter; it++) {
    s.interactionMetropolis();
    s.sampleTauInteractionPosterior();
    for(uword p = 0; p < P; p++) {
      tau2_trace(it, p) = s.tau2_interaction(p);
      double ss = 0.0;
      for(uword k = 0; k < K; k++) {
        for(uword b = 0; b < B; b++) {
          ss += s.gamma(p, k, b) * s.gamma(p, k, b);
        }
      }
      ss_trace(it, p) = ss;
    }
    for(uword k = 0; k < K; k++) {
      for(uword b = 0; b < B; b++) {
        gamma_trace(k, b, it) = s.gamma(0, k, b);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("tau2") = tau2_trace,
    Rcpp::Named("gamma_ss") = ss_trace,
    Rcpp::Named("gamma") = gamma_trace,
    Rcpp::Named("gamma_count") = conv_to<vec>::from(s.gamma_count)
  );
}
