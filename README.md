# batchmix
Semi-supervised and unsupervised Bayesian mixture models that 
simultaneously infer the cluster/class structure and a batch correction. 
Densities available are the multivariate normal and the multivariate t. 
The model sampler is implemented in C++. This package is aimed at analysis of 
low-dimensional data generated across several batches. See 
(Coleman et al. (2022))[https://doi.org/10.1101/2022.01.14.476352] for details
of the model.

The main functions a user should be aware of are ``runMCMCChains``, ``plotLikelihoods``, ``continueChains`` and ``processChains``. For an example of a workflow please see the short vignette.
