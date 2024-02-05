// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

arma::vec sample_cpp(const arma::vec &x,
                     const int size,
                     const bool replace,
                     Rcpp::NumericVector prob = Rcpp::NumericVector::create()) {
  arma::vec samp = Rcpp::RcppArmadillo::sample(x, size, replace, prob);
  return samp;
}

arma::uvec sample_cpp(const arma::uvec &x,
                      const int size,
                      const bool replace,
                      Rcpp::NumericVector prob = Rcpp::NumericVector::create()) {
  arma::uvec samp = Rcpp::RcppArmadillo::sample(x, size, replace, prob);
  return samp;
}
