% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSampledClusterMeans.R
\name{getSampledClusterMeans}
\alias{getSampledClusterMeans}
\title{Get sampled cluster means}
\usage{
getSampledClusterMeans(
  sampled_cluster_means,
  K = dim(sampled_cluster_means)[2],
  P = dim(sampled_cluster_means)[1],
  R = dim(sampled_cluster_means)[3],
  thin = 1
)
}
\arguments{
\item{sampled_cluster_means}{A 3D array of sampled cluster means.}

\item{K}{The number of clusters present. Defaults to the number of columns in
the batch mean matrix from the first sample.}

\item{P}{The dimension of the batch mean shifts. Defaults to the number of
rows in the batch mean matrix from the first sample.}

\item{R}{The number of iterations run. Defaults to the number of slices in
the sampled batch mean array.}

\item{thin}{The thinning factor of the sampler. Defaults to 1.}
}
\value{
A data.frame of three columns; the parameter, the sampled value and the iteration.
}
\description{
Given an array of sampled cluster means from the
``mixtureModel`` function, acquire a tidy version ready for ``ggplot2`` use.
}
\examples{

# Data in matrix format
X <- matrix(c(rnorm(100, 0, 1), rnorm(100, 3, 1)), ncol = 2, byrow = TRUE)

# Observed batches represented by integers
batch_vec <- sample(seq(1, 5), size = 100, replace = TRUE)

# MCMC iterations (this is too low for real use)
R <- 100
thin <- 5

# MCMC samples
samples <- runBatchMix(X, R, thin, batch_vec, "MVN")

batch_shift_df <- getSampledClusterMeans(samples$means, R = R, thin = thin)

}
