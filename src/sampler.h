// sampler.h
// =============================================================================
// include guard
#ifndef SAMPLER_H
#define SAMPLER_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>

// =============================================================================
// virtual sampler class


//' @name sampler
//' @title sampler
//' @description The virtual sampler class that is the parent to all specific 
//' implementations of the mixture model.
//' @field new Constructor \itemize{
//' \item Parameter: K - the number of components to model.
//' \item Parameter: B - the number of batches present.
//' \item Parameter: labels - N-vector of unsigned integers denoting initial 
//' clustering of the data .
//' \item Parameter: batch_vec - N-vector of unsigned integers denoting 
//' the observed grouping variable.
//' \item Parameter: concentration - K- vector of the prior hyperparameter for 
//' the class weights
//' \item Parameter: X - an N x P matrix of the observed data to model.
//' }
//' @field updateWeights Update the weights of each component based on current
//' clustering.
//' @field updateAllocation Sample a new clustering.
//' @field sampleFromPrior Sample parameter values from the prior distributions.
//' @field metropolisStep Perform Metropolis-Hastings sampling for the class and  
//' batch parameters.
//' @field calcBIC Calculate the BIC of the model.
//' @field logLikelihood Calculate the likelihood of a given data point in each
//' component. \itemize{
//' \item Parameter: x - a data point.
//' \item Parameter: b - the grouping variable category associated with x.
//' }
class sampler {

private:

public:


  arma::uword K = 0, B = 0, N = 0, P = 0, K_occ = 0, accepted = 0, N_fixed = 0;


  double observed_likelihood = 0.0,
    BIC = 0.0,
    complete_likelihood = 0.0;

  // fixed/unfixed_ind: which labels are observed (semi-supervised) vs. to be
  // sampled. A fixed vector of all zeroes gives the fully unsupervised model
  // - see updateAllocation().
  arma::uvec labels, N_k, batch_vec, N_b, KB_inds, B_inds, fixed, unfixed_ind;
  arma::vec concentration, ll, likelihood;
  arma::umat members;
  // alloc: the per-item allocation probability matrix. For fixed items this
  // holds a one-hot encoding of the known label; for unfixed items it is
  // populated by updateAllocation().
  arma::mat X, X_t, w, alloc;
  arma::field<arma::uvec> batch_ind;

  // Missing-data bookkeeping, shared by every concrete sampler (see
  // updateLatentData() in mvnSampler/mvtSampler/mvnSamplerSeparationStrategy,
  // and mvnSamplerMixed which additionally augments binary/censored entries
  // on top of this same NaN-driven baseline).
  //
  // X_raw/X_raw_t: the data exactly as supplied by the user - NaN at every
  // missing entry - and never modified after construction. X/X_t are the
  // *working* complete-data copy every likelihood/kernel actually reads;
  // a concrete sampler's updateLatentData() overwrites the missing entries
  // of X/X_t in place every MCMC sweep (see the *::updateLatentData()
  // implementations), but X_raw/X_raw_t are never touched again, and the
  // user's original R-level X object is never touched at all (Rcpp passes
  // this constructor a copy). Posterior draws of the imputed values are
  // exposed as a separate `latent_data` trace by the driver functions, not
  // by mutating X_raw.
  //
  // items_to_augment: indices of rows with at least one non-finite (NaN)
  // entry in X_raw, i.e. requiring imputation every sweep. Fully-observed
  // rows are skipped entirely by updateLatentData(), and if X has no
  // missing data at all this is empty, so updateLatentData() draws nothing
  // and every RNG call sequence - and hence every existing golden-master
  // test - is unaffected.
  arma::mat X_raw, X_raw_t;
  arma::uvec items_to_augment;

  // ===========================================================================
  // Batch x cluster interaction term in the mean (opt-in; see
  // sampleTauInteractionPosterior()/the *::interactionMetropolis()
  // overrides for the model).
  //
  // A first version of this put an independent N(0, tau2_interaction) prior
  // on every gamma(p, k, b) cell, relying purely on shrinkage (tight
  // tau2_interaction) for identifiability against the existing free mu_k/m_b
  // main effects. That does NOT work: mu_k + m_b + gamma_{k,b} is exactly
  // invariant under mu_k -> mu_k + c_k, m_b -> m_b + d_b, gamma_{k,b} ->
  // gamma_{k,b} - c_k - d_b for ANY c, d (the classical two-way ANOVA
  // main-effect/interaction confound, Gelman, 2005, "Analysis of variance -
  // why it is more important than ever", Annals of Statistics 33(1)), and
  // mu_k/m_b's own priors are far too diffuse to resist drifting along this
  // invariance - verified empirically: even an InvGamma(10, 0.1) hyperprior
  // on tau2_interaction (prior mean 0.011) left a single true nonzero
  // interaction cell unrecovered and every other cell sizeable, because
  // mu/m simply absorbed part of the signal that should have stayed additive.
  //
  // The fix used here is the standard one (effect coding / sum-to-zero
  // constrained interaction, e.g. Gelman 2005 Section 3; BDA3 Ch. 15):
  // gamma is confined for good to the subspace satisfying
  // sum_b gamma_{k,b} = 0 for every k and sum_k gamma_{k,b} = 0 for every b.
  // That subspace has exactly (K-1)(B-1) effective dimensions per feature,
  // and once gamma is restricted to it, mu_k + m_b is the unique best
  // additive explanation of the cell means and gamma is exactly (and only)
  // the part no additive mu_k + m_b could ever explain - the confound is
  // removed by construction, not by prior tightness. Both the prior
  // (sampleGammaPrior()) and every proposal (interactionMetropolis()) keep
  // gamma inside this subspace by construction (project a raw perturbation
  // via doubleCenterMatrix() before adding it), so it is never sampled
  // outside it and never needs correcting after the fact.
  bool include_interaction = false;
  double a_gamma = 2.0, b_gamma = 1.0, gamma_proposal_window = 0.1;
  arma::vec tau2_interaction;
  arma::cube gamma;
  arma::uvec gamma_count;

  // Projects an arbitrary K x B matrix onto the sum-to-zero subspace (every
  // row and every column sums to zero) via the standard two-way ANOVA
  // "double centering" identity M' = M - row_mean - col_mean + grand_mean.
  // Used to keep gamma's prior draws and proposals confined to that
  // subspace; see the section header comment above for why this matters.
  arma::mat doubleCenterMatrix(arma::mat M);

  // ===========================================================================
  // Batch-specific mixing weights (opt-in; see updateWeights()'s dispatch on
  // weight_prior_type). Rather than one global weight vector w shared by
  // every batch, each batch b can instead get its own weight vector
  // w_batch(b, .), built from the same additive-log-ratio (ALR)
  // parameterisation either way; what differs is the prior tying the K-1
  // free ALR coordinates together across batches:
  //
  //   weight_prior_type == 0 ("global", the default): no batch-specific
  //     weights at all - a single w shared by every batch, exactly the
  //     original behaviour. w_batch/eta_alr are allocated but never
  //     updated (they stay at their inert uniform/zero initial value).
  //
  //   weight_prior_type == 1 ("partial pooling" / exchangeable): each
  //     batch's ALR coordinate is a draw eta_{b,j} ~ N(mu_j, tau2_j) from a
  //     shared, ESTIMATED population mean mu_j and variance tau2_j, with no
  //     assumed order, distance or covariance structure between batches at
  //     all - batches are exchangeable (Gelman & Hill, 2007, "Data Analysis
  //     Using Regression and Multilevel/Hierarchical Models", ch. 12, for
  //     the general no-pooling/partial-pooling/complete-pooling framework
  //     this sits inside). tau2_j controls how much pooling happens: small
  //     tau2_j pulls every batch strongly towards the shared mu_j (close to
  //     "global"), large tau2_j lets batches vary close to independently.
  //     Unlike the interaction term, there is no other parameter competing
  //     to explain the same signal, so no identifiability fix is needed
  //     here beyond the ordinary conjugate hierarchical-model machinery.
  //
  //   weight_prior_type == 2 ("gp"): the ALR coordinates are linked across
  //     batches by a Gaussian process over an ordered/spatial
  //     batch_coordinates covariate (maternKernel32() +
  //     multinomialLogitGPLogKernel() in genericFunctions.h/.cpp), for
  //     batches with a genuine, known ordering in time or space where nearby
  //     batches are expected to be more similar than distant ones. Use
  //     "partial pooling" instead whenever that ordering/distance isn't
  //     known or isn't the point - it estimates only a common mean and
  //     variance, not a length-scale.
  //
  //     Each free ALR coordinate j has its own estimated intercept
  //     gp_beta(j) - eta_{.,j} ~ GP(gp_beta(j), K(batch_coordinates)) - not
  //     fixed at 0: Ren, Du, Carin & Dunson (2011), the paper this
  //     construction is otherwise built from, itself decomposes the
  //     logit as z(x)^T beta + w(x) with w ~ GP(0, K), i.e. an intercept
  //     plus a zero-mean deviation process, not a bare zero-mean GP - the
  //     latter would force every coordinate to revert to an equal-weight
  //     (1/K) split wherever the data are uninformative or far from every
  //     other batch, regardless of the batches' actual overall rate,
  //     unlike "partial pooling"'s pp_mu playing the equivalent role
  //     there. gp_beta(j) ~ N(0, pp_mu_prior_sd^2) a priori (the same
  //     weakly-informative scale already used for pp_mu, since the two
  //     play an analogous role) and is updated by a conjugate
  //     generalised-least-squares Gibbs step each sweep (see
  //     updateGPWeights()).
  //
  //     gp_length_scale's prior (gp_length_scale_prior_shape/rate, an
  //     Inverse-Gamma - NOT the log-normal an earlier version of this
  //     package used) is NOT the fixed constant it defaults to below once
  //     batch_coordinates has at least 2 distinct values:
  //     initialiseBatchWeightPrior() overwrites it with an empirical-Bayes,
  //     "boundary-avoiding" calibration to batch_coordinates' own scale via
  //     invGammaQuantileMatch() - see its implementation and
  //     genericFunctions.h for the full derivation and references
  //     (Betancourt, 2020, "Robust Gaussian Process Modeling," Stan case
  //     study; Fuglstad, Simpson, Lindgren & Rue, 2019, JASA, penalised-
  //     complexity priors for GP range/variance). A length-scale prior
  //     fixed in absolute units regardless of whether batch_coordinates
  //     are e.g. small integer indices or real-valued collection times
  //     spanning months/years is a scale mismatch waiting to happen - and
  //     a log-normal (rather than Inverse-Gamma) calibration, tried in an
  //     earlier version of this fix, turned out to be its own trap: one
  //     coincidentally-close pair of batches among otherwise-spread-out
  //     ones inflates the log-normal's SD enough that its heavy right tail
  //     puts real prior mass on length scales tens to hundreds of times the
  //     whole batch_coordinates range, and - because the likelihood is
  //     essentially flat once length-scale already exceeds that range - the
  //     sampler can wander arbitrarily far up that flat ridge and get stuck
  //     there, collapsing the entire GP to one shared value regardless of
  //     any real per-batch structure. Inverse-Gamma's much faster-decaying
  //     right tail (the reason it, not log-normal, is the standard
  //     recommendation for this specific prior) keeps that from happening.
  //
  //     The GP prior's quadratic form is evaluated via a triangular solve
  //     against gp_chol (the covariance's Cholesky factor), never via an
  //     explicit matrix inverse - see multinomialLogitGPLogKernel() in
  //     genericFunctions.h for why this, not gp_cov's inverse, is the
  //     numerically-preferred modern approach; gp_cov_inv does not exist as
  //     a field for this reason. Similarly, gpHyperparameterMetropolis()
  //     proposes new (tau2, length_scale) via a non-centred/whitened
  //     reparameterisation of eta_alr specifically for that step (see its
  //     own documentation) - the standard fix for the otherwise-severe
  //     "funnel" coupling between a hierarchical model's variance/
  //     length-scale and its latent values when both are updated in their
  //     natural, entangled ("centred") form (Neal, 2003; Papaspiliopoulos,
  //     Roberts & Skold, 2007; Betancourt & Girolami, 2015; Betancourt,
  //     2020, applying the same fix to a GP specifically) - rather than the
  //     plain, funnel-prone joint update an earlier version of this fix
  //     used.
  arma::uword weight_prior_type = 0;
  bool sample_gp_hyperparameters = false;
  double gp_tau2 = 1.0, gp_length_scale = 1.0, gp_jitter = 1e-6,
    eta_proposal_window = 0.1, gp_hyperparameter_proposal_window = 0.1,
    gp_tau2_prior_shape = 2.0, gp_tau2_prior_rate = 4.0,
    gp_length_scale_prior_shape = 2.0, gp_length_scale_prior_rate = 2.0,
    pp_tau2_prior_shape = 2.0, pp_tau2_prior_rate = 1.0,
    pp_mu_prior_sd = 10.0;
  arma::vec batch_coordinates, pp_mu, pp_tau2, gp_beta;
  arma::mat gp_cov, gp_chol, w_batch, eta_alr;
  arma::uvec eta_count;
  arma::uword gp_hyperparameter_count = 0;

  // Parametrised class
  sampler(
    arma::uword _K,
    arma::uword _B,
    arma::uvec _labels,
    arma::uvec _batch_vec,
    arma::vec _concentration,
    arma::mat _X,
    arma::uvec _fixed);

  // Destructor (I haven't found a way of defining destructors outside headers)
  virtual ~sampler() { };

  // Functions required of all mixture models
  // Generic functions
  virtual void updateWeights();
  virtual void updateAllocation();

  // Interaction-term setup/priors (shared across every concrete sampler,
  // since they only touch K/B/P/gamma/tau2_interaction, all declared here).
  // Called from each subclass's constructor body once K/B/P are known.
  void initialiseInteraction(bool _include_interaction, double _gamma_proposal_window,
                              double _a_gamma, double _b_gamma);
  void sampleGammaPrior();
  void sampleTauInteractionPosterior();

  // Batch-weight-prior setup/updates (shared: only touches
  // K/B/N_b/labels/batch_vec/w/gp_*/pp_*/eta_alr/w_batch, all declared here).
  void initialiseBatchWeightPrior(arma::uword _weight_prior_type,
                                   arma::vec _batch_coordinates,
                                   double _gp_tau2,
                                   double _gp_length_scale,
                                   double _eta_proposal_window,
                                   bool _sample_gp_hyperparameters,
                                   double _gp_hyperparameter_proposal_window,
                                   double _pp_tau2_prior_shape,
                                   double _pp_tau2_prior_rate,
                                   double _pp_mu_prior_sd);
  arma::mat computeBatchClassCounts();
  void updateSimplexFromALR(arma::uword n_free);
  void updateGPWeights();
  void updatePartialPoolingWeights();
  void gpHyperparameterMetropolis();

  // Builds a Matern-3/2 GP covariance matrix and its LOWER-triangular
  // Cholesky factor for a candidate (tau2, length_scale) pair, escalating
  // the diagonal jitter geometrically from a tau2-scaled starting point
  // until the Cholesky factorisation succeeds (Armadillo's chol() returns
  // false rather than throwing on failure here) - the standard modern
  // safeguard for covariance matrices that can be arbitrarily
  // ill-conditioned depending on the input locations (GPyTorch's/GPflow's
  // escalating default jitter, Gardner, Pleiss, Bindel, Weinberger &
  // Wilson, 2018, "GPyTorch," NeurIPS), used because no fixed proportional
  // jitter is safe for every possible batch_coordinates configuration (see
  // initialiseBatchWeightPrior()). Deliberately factorises via chol(), not
  // inv_sympd(): every consumer of the result (multinomialLogitGPLogKernel(),
  // gpHyperparameterMetropolis()'s whitening step) only ever needs a
  // triangular solve against the Cholesky factor, never the covariance's
  // explicit inverse - see genericFunctions.h for why that is the
  // numerically-preferred modern approach. Shared by
  // initialiseBatchWeightPrior() and gpHyperparameterMetropolis(), the two
  // places a covariance is built from scratch for a candidate (tau2,
  // length_scale).
  void buildWellConditionedGPChol(
    double tau2, double length_scale,
    arma::mat& cov, arma::mat& chol_factor, double& jitter_used
  );

  // Extra free-parameter count (beyond each concrete calcBIC()'s baseline
  // cluster/batch terms) contributed by the opt-in interaction term and/or
  // batch-specific weight prior, when either is switched on. Shared across
  // every concrete sampler since it only touches members declared here
  // (K_occ/B/P/include_interaction/weight_prior_type/
  // sample_gp_hyperparameters) - see calcBIC() in each subclass, and
  // sampleTauInteractionPosterior() above for the (K-1)(B-1)-per-feature
  // effective dimension of gamma this mirrors.
  double structuralExtraBICParams() const;

  // Mixture specific functions (therefore virtual at this level)
  virtual void metropolisStep() = 0;
  virtual void sampleFromPriors() = 0;
  virtual void calcBIC() = 0;
  virtual arma::vec itemLogLikelihood(arma::vec x, arma::uword b) = 0;

};

#endif /* SAMPLER_H */