# BatchMixtureModel
Bayesian mixture models that explicitly infer both the group/class/cluster parameters and the batch parameters. These are used to define a partition of the batch-corrected dataset.

The main functions a user should be aware of are ``runMCMCChains``, ``plotLikelihoods``, ``continueChains`` and ``processChains``. For an example of a workflow please see the short vignette.
