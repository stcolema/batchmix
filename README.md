# batchmix

Downloads from CRAN in the past month:

[![](https://cranlogs.r-pkg.org/badges/batchmix)](https://cran.r-project.org/package=batchmix)

Semi-supervised and unsupervised Bayesian mixture models that 
simultaneously infer the cluster/class structure and a batch correction. 
Densities available are the multivariate normal and the multivariate t. 
The model sampler is implemented in C++. This package is aimed at analysis of 
low-dimensional data generated across several batches. See 
[Coleman et al. (2022)](https://doi.org/10.1101/2022.01.14.476352) for details
of the model.

## Advice on using the package

The main functions a user should be aware of are ``runMCMCChains``, 
``plotLikelihoods``, ``plotAcceptanceRates``, ``continueChains`` and 
``processChains``. 

Parameters are sampled using Metropolis-Hastings so checking that the 
acceptance rate is important. We recommend aiming for acceptance
rates between 0.1 and 0.5 for the class and batch means and batch 
scales ($\mu_k$, $m_b$ and $S_b$ respectively). In our testing, an 
acceptance rate of at least 0.4 for the class covariance matrices
tended to suggest the sampler is exploring well, but smaller values 
were frequently associated with poor behaviour. The degrees of freedom tend 
to have very high acceptance rates in our testing regardless of the 
sampling window.

We recommend running a small number of chains for a small number of 
iterations to assess the acceptance rates before committing the 
computational resourcces to run a full analysis.

For an example of a workflow please see the short vignette.
