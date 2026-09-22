// mvnSamplerMixed.h
// =============================================================================
// include guard
#ifndef MVNSAMPLERMIXED_H
#define MVNSAMPLERMIXED_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "sampler.h"
# include "genericFunctions.h"
# include "mvnSamplerSeparationStrategy.h"

// =============================================================================
// mvnSamplerMixed class
//
// Extends the LKJ-based separation-strategy sampler (mu, m batch shift,
// R/sigma covariance decomposition, S batch scale) with a shared latent
// Gaussian layer that supports:
//
//   - continuous columns, exactly as in mvnSamplerSeparationStrategy;
//   - binary columns observed through a probit link (Albert & Chib, 1993,
//     JASA; Chib & Greenberg, 1998, Biometrika, for the multivariate case);
//   - missing-at-random entries in continuous columns;
//   - left- or right-censored entries in continuous columns (e.g. assay
//     values below/above a detection limit), recorded at the known
//     censoring bound rather than as missing.
//
// This is the "one shared latent Gaussian, different observation models
// per column" construction of Dunson (2000, JRSS-B) and Dunson & Herring
// (2005, Biostatistics): every column is generated from the same
// underlying (mu_k + m_b, Sigma_k) latent structure, and the augmentation
// step below reconstructs a complete-data draw of that latent vector from
// whatever is actually observed in each column, before every other
// parameter update proceeds exactly as in the fully-continuous sampler.
//
// Identifiability: a multivariate probit model is only identified up to
// the CORRELATION matrix of the latent vector - rescaling a coordinate's
// latent variance leaves sign(z) (hence the observed binary outcome)
// unchanged (Chib & Greenberg, 1998). Both scale layers of the R/sigma
// decomposition - the cluster-level marginal sd (sigma) and the batch-level
// multiplicative scale (S) - are therefore fixed at 1 for binary columns
// and never touched by their respective Metropolis/prior-sampling steps;
// only the shared correlation matrix R and the location parameters
// (mu, m) are estimated for those columns. Continuous columns are
// unaffected and keep their usual free sigma/S.
class mvnSamplerMixed: public mvnSamplerSeparationStrategy {

public:

  // 0 = continuous column, 1 = binary column observed via a probit link.
  arma::uvec column_type;

  // Per-entry censoring code for continuous columns: 0 = not censored,
  // 1 = left-censored (true value < X_raw, the recorded value being the
  // known lower detection/quantification bound), 2 = right-censored (true
  // value > X_raw, recorded value the known upper bound). Ignored for
  // binary columns and for missing (NaN) entries.
  arma::umat censor_code;

  // The data as supplied by the user: NaN for missing entries, the
  // recorded bound for censored entries, {0, 1} for binary columns,
  // otherwise the observed value. Unlike the inherited X/X_t (which this
  // class overwrites every sweep with the current complete-data draw of
  // the shared latent vector), X_raw/X_raw_t never change after
  // construction.
  arma::mat X_raw, X_raw_t;

  // Indices of items requiring at least one entry to be augmented (i.e.
  // that have a missing, censored, or binary entry); items that are fully
  // observed and continuous are skipped entirely in updateLatentData().
  arma::uvec items_to_augment;

  // A simple column-mean imputation of X, used ONLY to give the base
  // class's empirical-Bayes prior-hyperparameter calculations (mu_0,
  // global_cov, etc., which use mean()/cov() and do not skip NaN) a
  // finite matrix to work with; it must run before the base class
  // constructors do, i.e. inside this class's member-initializer list,
  // which is why it is static (no `this` exists yet at that point). It
  // has no bearing on the actual data augmentation carried out by
  // updateLatentData(), which is driven entirely by X_raw/column_type/
  // censor_code.
  static arma::mat imputeForPriorSetup(arma::mat X);

  mvnSamplerMixed(
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
  );

  // Destructor
  virtual ~mvnSamplerMixed() { };

  // The Gibbs step that reconstructs a complete-data draw z_i for every
  // item requiring augmentation, given the current parameter values, and
  // writes the result into the inherited X/X_t fields. Must be called
  // once per MCMC sweep, before updateAllocation()/metropolisStep(), so
  // that every other (unmodified, inherited) step operates on complete
  // continuous-scale data exactly as in mvnSamplerSeparationStrategy.
  void updateLatentData();

  // Overridden to hold sigma_{k,p} = 1 and S_{b,p} = 1 fixed for binary
  // columns p (see class-level identifiability note above).
  void sampleCovPrior();
  void sampleSPrior();
  void sigmaMHStep();
  void batchScaleMetropolis();

  // Overridden: binary/probit columns have sigma and S fixed at 1 rather
  // than estimated, so they contribute fewer free parameters per cluster/
  // batch than a continuous column does. See the .cpp for the derivation,
  // and an important caveat about what observed_likelihood represents
  // (the augmented-data likelihood, not the true marginal p(y|theta))
  // whenever any column is binary.
  void calcBIC();

};

#endif /* MVNSAMPLERMIXED_H */
