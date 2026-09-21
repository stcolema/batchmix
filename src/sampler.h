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


  arma::uword K = 0, B = 0, N = 0, P = 0, K_occ = 0, accepted = 0;


  double observed_likelihood = 0.0,
    BIC = 0.0,
    complete_likelihood = 0.0;

  arma::uvec labels, N_k, batch_vec, N_b, KB_inds, B_inds;
  arma::vec concentration, ll, likelihood;
  arma::umat members;
  arma::mat X, X_t, w; //, alloc;
  arma::field<arma::uvec> batch_ind;

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
  //   weight_prior_type == 2 ("gp"): as before - the ALR coordinates are
  //     linked across batches by a Gaussian process over an ordered/spatial
  //     batch_coordinates covariate (squaredExponentialKernel() +
  //     multinomialLogitGPLogKernel() in genericFunctions.h/.cpp), for
  //     batches with a genuine, known ordering in time or space where nearby
  //     batches are expected to be more similar than distant ones. Use
  //     "partial pooling" instead whenever that ordering/distance isn't
  //     known or isn't the point - it estimates only a common mean and
  //     variance, not a length-scale.
  arma::uword weight_prior_type = 0;
  bool sample_gp_hyperparameters = false;
  double gp_tau2 = 1.0, gp_length_scale = 1.0, gp_jitter = 1e-6,
    eta_proposal_window = 0.1, gp_hyperparameter_proposal_window = 0.1,
    gp_tau2_prior_shape = 2.0, gp_tau2_prior_rate = 1.0,
    gp_length_scale_prior_mean = 0.0, gp_length_scale_prior_sd = 1.0,
    pp_tau2_prior_shape = 2.0, pp_tau2_prior_rate = 1.0,
    pp_mu_prior_sd = 10.0;
  arma::vec batch_coordinates, pp_mu, pp_tau2;
  arma::mat gp_cov, gp_cov_inv, w_batch, eta_alr;
  arma::uvec eta_count;
  arma::uword gp_hyperparameter_count = 0;

  // Parametrised class
  sampler(
    arma::uword _K,
    arma::uword _B,
    arma::uvec _labels,
    arma::uvec _batch_vec,
    arma::vec _concentration,
    arma::mat _X);

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
  double gpHyperparameterLogKernel(arma::mat cov, arma::mat cov_inv);
  void gpHyperparameterMetropolis();

  // Mixture specific functions (therefore virtual at this level)
  virtual void metropolisStep() = 0;
  virtual void sampleFromPriors() = 0;
  virtual void calcBIC() = 0;
  virtual arma::vec itemLogLikelihood(arma::vec x, arma::uword b) = 0;

};

#endif /* SAMPLER_H */