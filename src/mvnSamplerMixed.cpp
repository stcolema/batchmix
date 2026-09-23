// mvnSamplerMixed.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "pdfs.h"
# include "sampler.h"
# include "mvnSamplerMixed.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// mvnSamplerMixed class

// See the header: this exists purely to give the base class's NaN-
// intolerant empirical-Bayes prior calculations (mean(), cov() over the
// raw data) a finite matrix, before this class's own constructor body has
// a chance to run and set up the real (properly augmented) starting data.
arma::mat mvnSamplerMixed::imputeForPriorSetup(arma::mat X) {
  return imputeColumnMeans(X);
};

mvnSamplerMixed::mvnSamplerMixed(
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
  double _eta,
  arma::uvec _column_type,
  arma::umat _censor_code
) : mvnSamplerSeparationStrategy(_K,
_B,
_mu_proposal_window,
_r_proposal_window,
_sigma_proposal_window,
_m_proposal_window,
_S_proposal_window,
_labels,
_batch_vec,
_concentration,
mvnSamplerMixed::imputeForPriorSetup(_X),
_fixed,
_m_scale,
_rho,
_theta,
_sample_m_scale,
_eta
) {

  column_type = _column_type;
  censor_code = _censor_code;
  X_raw = _X;
  X_raw_t = _X.t();

  // Initialise the working complete-data copy (the inherited X/X_t, which
  // this class overwrites every sweep): observed continuous values and
  // {0,1} binary values are kept as-is for now (sign-mapped below for
  // binary columns), censored entries start at their recorded bound, and
  // missing entries start at the observed column mean (continuous) or 0
  // (binary). updateLatentData() will move all of these to a proper draw
  // before they are ever used in a likelihood evaluation.
  X = X_raw;
  for(uword p = 0; p < P; p++) {
    if(column_type(p) == 1) {
      for(uword n = 0; n < N; n++) {
        if(std::isnan(X_raw(n, p))) {
          X(n, p) = 0.0;
        } else {
          X(n, p) = (X_raw(n, p) > 0.5) ? 1.0 : -1.0;
        }
      }
    } else {
      vec col_p = X_raw.col(p);
      uvec observed = find_finite(col_p);
      double col_mean = observed.n_elem > 0 ? mean(col_p.elem(observed)) : 0.0;
      for(uword n = 0; n < N; n++) {
        if(std::isnan(X_raw(n, p))) {
          X(n, p) = col_mean;
        }
      }
    }
  }
  X_t = X.t();

  // Items with nothing to augment (fully observed, continuous, not
  // censored) are skipped entirely by updateLatentData().
  std::vector<uword> to_augment;
  bool needs = false;
  for(uword n = 0; n < N; n++) {
    needs = false;
    for(uword p = 0; p < P; p++) {
      if(column_type(p) == 1 || std::isnan(X_raw(n, p)) || censor_code(n, p) != 0) {
        needs = true;
        break;
      }
    }
    if(needs) {
      to_augment.push_back(n);
    }
  }
  items_to_augment = arma::conv_to<uvec>::from(to_augment);
};

// The shared-latent-Gaussian data augmentation step (Albert & Chib, 1993;
// Chib & Greenberg, 1998; Dunson, 2000). For every item with at least one
// missing, censored or binary entry, cycle through those entries in a
// systematic-scan Gibbs sweep: the conditional distribution of z_p given
// every other current entry of the same item is Normal, and is obtained
// cheaply from the already-computed precision matrix cov_comb_inv (no
// per-item matrix inversion needed), then either used directly (missing),
// or truncated to the half-line/interval implied by the binary outcome or
// the censoring bound.
void mvnSamplerMixed::updateLatentData() {

  vec z_i(P), eta(P);
  double cond_mean = 0.0, cond_var = 0.0, cond_sd = 0.0, lambda_pp = 0.0;
  uword k = 0, b = 0, kb = 0, code = 0;
  bool is_missing = false, is_binary = false;

  for(auto& n : items_to_augment) {

    k = labels(n);
    b = batch_vec(n);
    kb = k * B + b;

    z_i = X_t.col(n);
    eta = mean_sum.col(kb);

    for(uword p = 0; p < P; p++) {

      is_binary = (column_type(p) == 1);
      is_missing = std::isnan(X_raw_t(p, n));
      code = is_binary ? 0 : censor_code(n, p);

      if(!is_binary && !is_missing && code == 0) {
        continue;
      }

      lambda_pp = cov_comb_inv(p, p, kb);
      cond_var = 1.0 / lambda_pp;
      cond_mean = eta(p) - cond_var * (
        arma::dot(cov_comb_inv.slice(kb).row(p), z_i - eta) - lambda_pp * (z_i(p) - eta(p))
      );
      cond_sd = std::sqrt(cond_var);

      if(is_binary) {
        if(is_missing) {
          z_i(p) = cond_mean + cond_sd * randn();
        } else if(X_raw_t(p, n) > 0.5) {
          z_i(p) = rTruncNorm(cond_mean, cond_sd, 0.0, arma::datum::inf);
        } else {
          z_i(p) = rTruncNorm(cond_mean, cond_sd, -arma::datum::inf, 0.0);
        }
      } else if(is_missing) {
        z_i(p) = cond_mean + cond_sd * randn();
      } else if(code == 1) {
        // Left-censored: the true value is below the recorded bound.
        z_i(p) = rTruncNorm(cond_mean, cond_sd, -arma::datum::inf, X_raw_t(p, n));
      } else if(code == 2) {
        // Right-censored: the true value is above the recorded bound.
        z_i(p) = rTruncNorm(cond_mean, cond_sd, X_raw_t(p, n), arma::datum::inf);
      }
    }

    X_t.col(n) = z_i;
  }

  X = X_t.t();
};

// Fixes sigma_{k,p} = 1 for binary columns (see class-level identifiability
// note in the header); continuous columns are sampled exactly as in the
// parent class.
void mvnSamplerMixed::sampleCovPrior() {
  vec log_sigma(P);
  log_sigma.zeros();
  for(uword k = 0; k < K; k++){
    R.slice(k) = sampleLKJCorrelationMatrix(P, eta);
    log_sigma = randn<vec>(P, distr_param(beta, xi));

    sigma.col(k) = exp(log_sigma);
    for(uword p = 0; p < P; p++) {
      if(column_type(p) == 1) {
        sigma(p, k) = 1.0;
      }
    }
    Sigma_mat.slice(k).diag() = sigma.col(k);

    cov.slice(k) = Sigma_mat.slice(k) * R.slice(k) * Sigma_mat.slice(k);
  }
};

// Fixes S_{b,p} = 1 for binary columns.
void mvnSamplerMixed::sampleSPrior() {
  for(uword b = 0; b < B; b++){
    for(uword p = 0; p < P; p++){
      if(column_type(p) == 1) {
        S(p, b) = 1.0;
      } else {
        S(p, b) = S_loc + 1.0 / randg<double>( distr_param(rho, 1.0 / theta) );
      }
    }
  }
};

// As mvnSamplerSeparationStrategy::batchScaleMetropolis(), but S_{b,p} is
// never proposed (and so never moves from 1) for binary columns p.
void mvnSamplerMixed::batchScaleMetropolis() {

  bool next = false;
  double u = 0.0, proposed_model_score = 0.0, acceptance_prob = 0.0, current_model_score = 0.0;
  vec S_proposed(P), proposed_cov_comb_log_det(K);
  cube proposed_cov_comb(P, P, K), proposed_cov_comb_inv(P, P, K);

  proposed_cov_comb_log_det.zeros();
  proposed_cov_comb.zeros();
  proposed_cov_comb_inv.zeros();

  for(uword b = 0; b < B; b++) {

    next = false;
    acceptance_prob = 0.0, proposed_model_score = 0.0, current_model_score = 0.0;
    proposed_cov_comb.zeros();
    S_proposed = S.col(b);

    for(uword p = 0; p < P; p++) {
      if(column_type(p) == 1) {
        continue;
      }
      if((S(p, b) - S_loc) * S_proposal_window < 0.0){
        Rcpp::stop("\n\nCurent batch scale equals S_loc");
      }
      S_proposed(p) = S_loc + randg( distr_param( (S(p, b) - S_loc) * S_proposal_window, 1.0 / S_proposal_window) );

      if(S_proposed(p) <= S_loc) {
        next = true;
      }

      // See mvnSampler::batchScaleMetropolis for the derivation/empirical
      // check: reverse density q(current|proposed) -> proposed_model_score,
      // forward density q(proposed|current) -> current_model_score.
      proposed_model_score += gammaLogLikelihood(S(p, b) - S_loc, (S_proposed(p) - S_loc) * S_proposal_window, S_proposal_window);
      current_model_score += gammaLogLikelihood(S_proposed(p) - S_loc, (S(p, b) - S_loc) * S_proposal_window, S_proposal_window);
    }

    if(next) {
      continue;
    }

    proposed_cov_comb = cov;
    // See mvnSamplerSeparationStrategy::batchScaleMetropolis() for why this
    // is guarded rather than assumed: cov.slice(k) can be PD in exact
    // arithmetic yet numerically singular in floating point, and
    // S_proposed only inflates the diagonal, so it cannot fix that - treat
    // a numerically-degenerate proposal as an automatic reject.
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

    proposed_model_score += sLogKernel(b, S_proposed, proposed_cov_comb_log_det, proposed_cov_comb_inv);
    current_model_score += sLogKernel(b, S.col(b), cov_comb_log_det.col(b), cov_comb_inv.slices(KB_inds + b));

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

// As mvnSamplerSeparationStrategy::sigmaMHStep(), but sigma_{k,p} is never
// proposed (and so never moves from 1) for binary columns p.
void mvnSamplerMixed::sigmaMHStep() {
  bool next = false;

  double u = 0.0,
    proposed_model_score = 0.0,
    acceptance_prob = 0.0,
    current_model_score = 0.0;

  arma::vec sigma_proposed(P), proposed_cov_comb_log_det(B);
  arma::mat sigma_mat_proposed(P, P), proposed_cov(P, P), proposed_cov_inv(P, P);
  arma::cube proposed_cov_comb(P, P, B), proposed_cov_comb_inv(P, P, B);

  proposed_cov.zeros();
  proposed_cov_inv.zeros();
  proposed_cov_comb_log_det.zeros();
  proposed_cov_comb.zeros();
  proposed_cov_comb_inv.zeros();

  for(arma::uword k = 0; k < K ; k++) {

    sigma_mat_proposed.diag() = sigma.col(k);
    sigma_proposed = sigma.col(k);
    next = false;

    proposed_cov_comb.zeros();

    acceptance_prob = 0.0, proposed_model_score = 0.0, current_model_score = 0.0;

    if(N_k(k) > 0) {
      for(uword p = 0; p < P; p++) {

        if(column_type(p) == 1) {
          continue;
        }

        sigma_proposed(p) = randg( distr_param( sigma(p, k) * sigma_proposal_window, 1.0 / sigma_proposal_window) );

        if(sigma_proposed(p) <= 1e-8) {
          next = true;
        }

        // Cross-referenced asymmetric proposal density (see
        // mvnSamplerSeparationStrategy::sigmaMHStep for the derivation) -
        // NOT each value scored under its own proposal shape. Reverse density
        // q(current|proposed) -> proposed_model_score, forward density
        // q(proposed|current) -> current_model_score.
        proposed_model_score += gammaLogLikelihood(sigma(p, k), sigma_proposed(p) * sigma_proposal_window, sigma_proposal_window);
        current_model_score += gammaLogLikelihood(sigma_proposed(p), sigma(p, k) * sigma_proposal_window, sigma_proposal_window);
      }

      if(next) {
        continue;
      }
    }

    // sigma_proposed/R.slice(k) can each be individually valid (fresh
    // LogNormal draw / RW-Gamma proposal; valid correlation matrix) yet
    // their product PD in exact arithmetic but numerically singular in
    // floating point (confirmed empirically for
    // mvnSamplerSeparationStrategy::sigmaMHStep()/rMHStep(), which this
    // mirrors). For N_k(k) == 0 (an unconditional prior draw, not an MH
    // proposal to accept/reject), redraw sigma_proposed - a
    // negligible-measure rejection of the pathological tail, not a bias on
    // the prior - up to max_attempts times; for N_k(k) > 0 (a genuine MH
    // proposal), a single numerically-degenerate outcome is simply an
    // automatic reject, matching `next`/`continue` above.
    bool cov_ok = false;
    const uword max_attempts = (N_k(k) == 0) ? 50 : 1;
    for(uword attempt = 0; attempt < max_attempts && !cov_ok; attempt++) {

      if(N_k(k) == 0) {
        // Sample fresh from the prior for continuous columns (matching
        // sampleCovPrior()'s log_sigma ~ N(beta, xi), sigma = exp(log_sigma)),
        // rather than an uncorrected, force-accepted random walk off the
        // cluster's current (possibly stale) sigma - see
        // mvnSamplerSeparationStrategy::sigmaMHStep() for the same fix.
        // Binary columns' sigma stays fixed at 1 for identifiability.
        for(uword p = 0; p < P; p++) {
          if(column_type(p) == 1) {
            continue;
          }
          sigma_proposed(p) = std::exp(arma::randn(distr_param(beta, xi)));
        }
      }

      sigma_mat_proposed.diag() = sigma_proposed;
      proposed_cov = sigma_mat_proposed * R.slice(k) * sigma_mat_proposed;
      cov_ok = arma::inv_sympd(proposed_cov_inv, proposed_cov);
      if(!cov_ok) {
        continue;
      }

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
    }
    if(!cov_ok) {
      continue;
    }

    if(N_k(k) > 0) {
      proposed_model_score += sigmaLogKernel(k,
                                             proposed_cov_comb_log_det,
                                             sigma_proposed,
                                             proposed_cov,
                                             proposed_cov_inv,
                                             proposed_cov_comb_inv
      );

      current_model_score += sigmaLogKernel(k,
                                        cov_comb_log_det.row(k).t(),
                                        sigma.col(k),
                                        cov.slice(k),
                                        cov_inv.slice(k),
                                        cov_comb_inv.slices(k * B + B_inds)
      );

      u = arma::randu();
      acceptance_prob = std::min(1.0, std::exp(proposed_model_score - current_model_score));
    }

    if( (u < acceptance_prob) || (N_k(k) == 0) ){
      sigma_count(k)++;
      sigma.col(k) = sigma_proposed;
      Sigma_mat.slice(k).diag() = sigma_proposed;

      cov.slice(k) = proposed_cov;
      cov_inv.slice(k) = proposed_cov_inv;
      cov_log_det(k) = arma::log_det(proposed_cov).real();

      for(arma::uword b = 0; b < B; b++) {
        cov_comb.slice(k * B + b) = proposed_cov_comb.slice(b);
        cov_comb_log_det(k, b) = proposed_cov_comb_log_det(b);
        cov_comb_inv.slice(k * B + b) = proposed_cov_comb_inv.slice(b);
      }
    }
  }
};

// Per cluster: 1 (weight) + P (mean; free for every column, including
// binary ones - only the SCALE of a binary column's latent variance is
// unidentified, not its location) + P*(P-1)/2 (R's free off-diagonal
// entries; also free for every column pair, including binary-binary and
// binary-continuous pairs - correlations involving a probit-linked
// column are identified, this is the standard tetrachoric/polyserial
// correlation result) + n_continuous (sigma; only continuous columns
// have a free marginal SD, binary columns are fixed at 1).
//
// Per batch: P (shift; free for every column) + n_continuous (scale;
// only continuous columns are free).
//
// IMPORTANT CAVEAT: unlike the fully-continuous samplers, observed_likelihood
// here is the log-density of the augmented complete data (X_t, including
// the latent draws for binary/missing/censored entries), not the true
// marginal log p(y | theta) of the observed outcomes - that would require
// integrating out the augmented latent variables, e.g. via GHK simulation
// of the relevant multivariate normal orthant probabilities, which is not
// implemented. BIC values from this sampler are therefore only a rough,
// within-model convergence/complexity diagnostic; they should not be
// used to compare against BIC from a model with no binary/missing/
// censored columns, or trusted for rigorous model selection, without
// that correction.
void mvnSamplerMixed::calcBIC() {

  double n_continuous = (double) accu(column_type == 0);
  double n_param_cluster_mixed = 1.0 + P + P * (P - 1) * 0.5 + n_continuous;
  double n_param_batch_mixed = P + n_continuous;

  // structuralExtraBICParams() adds the interaction term's and/or the
  // batch-specific weight prior's extra parameters when either is enabled
  // (0 in the default configuration) - see sampler.h/.cpp.
  BIC = 2 * observed_likelihood - (n_param_cluster_mixed * K_occ + n_param_batch_mixed * B + structuralExtraBICParams()) * std::log(N);

};
