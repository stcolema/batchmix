// pdfs.h
// =============================================================================
// include guard
#ifndef PDFS_H
#define PDFS_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>

// =============================================================================
// Log-likelihood functions used in Metropolis-hastings

//' @title Gamma log-likelihood
//' @description Used in calculating model probability in Metropolis-Hastings 
//' algorithm when proposals are from the Gamma distribution.
//' @param x - double; the value to calculate the unnormalised likelihood of.
//' @param shape - double; the shape of the Gamma distribution.
//' @param rate - double; the rate of the Gamma distribution
//' @return the unnormalised log-likelihood of x in a Gamma with parameters shape 
//' and rate.
// [[Rcpp::export]]
double gammaLogLikelihood(double x, double shape, double rate);

//' @title Inverse gamma log-likelihood
//' @description Used in calculating model probability in Metropolis-Hastings 
//' algorithm when proposals are from the inverse-Gamma distribution.
//' @param x - double; the value to calculate the likelihood of.
//' @param shape - double; the shape of the inverse-Gamma distribution.
//' @param scale - double; the scale of the inverse-Gamma distribution
//' @return the unnormalised log-likelihood of x in a inverse-Gamma with parameters 
//' shape and scale.
// [[Rcpp::export]]
double invGammaLogLikelihood(double x, double shape, double scale);

//' @title Wishart log-likelihood
//' @description Used in calculating model probability in Metropolis-Hastings 
//' algorithm when proposals are from the Wishart distribution.
//' @param X - matrix; the matrix to calculate the likelihood of.
//' @param V - matrix; the scale of the Wishart distribution.
//' @param n - double; the degrees of freedom for the Wishart distribution.
//' @param P - unsigned integer; the dimension of X.
//' @return the unnormalised log-likelihood of X in a Wishart with parameters V 
//' and n.
// [[Rcpp::export]]
double wishartLogLikelihood(arma::mat X, arma::mat V, double n, arma::uword P);

//' @title Inverse-Wishart log-likelihood
//' @description Used in calculating model probability in Metropolis-Hastings 
//' algorithm when proposals are from the Wishart distribution.
//' @param X - matrix; the matrix to calculate the likelihood of.
//' @param Psi - matrix; the scale of the inverse-Wishart distribution.
//' @param nu - double; the degrees of freedom for the inverse-Wishart distribution.
//' @param P - unsigned integer; the dimension of X.
//' @return the unnormalised log-likelihood of X in a inverse-Wishart with parameters Psi 
//' and nu.
// [[Rcpp::export]]
double invWishartLogLikelihood(arma::mat X, arma::mat Psi, double nu, arma::uword P);

//' @title LKJ correlation matrix log-likelihood
//' @description The unnormalised log-density of the LKJ(eta) distribution
//' for a correlation matrix R (Lewandowski, Kurowicka & Joe, 2009,
//' J. Multivariate Analysis 100(9)). The normalising constant depends only
//' on eta and the dimension P, not on R, so it is safe to omit whenever eta
//' is held fixed (e.g. in a Metropolis-Hastings acceptance ratio); it is
//' dropped here.
//' @param R - matrix; a P x P correlation matrix (symmetric, unit diagonal,
//' positive definite).
//' @param eta - double; the LKJ concentration/shape parameter. eta = 1 gives
//' a uniform density over the space of correlation matrices; eta > 1
//' concentrates mass towards the identity (weaker correlations); eta < 1
//' concentrates mass away from the identity (stronger correlations).
//' @return the unnormalised log-density of R under LKJ(eta).
// [[Rcpp::export]]
double lkjLogLikelihood(arma::mat R, double eta);

#endif /* PDFS_H */