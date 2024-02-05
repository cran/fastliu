#ifndef samp_H
#define samp_H

#include <RcppArmadillo.h>
arma::vec sample_cpp(const arma::vec &x,
           const int size,
           const bool replace,
           Rcpp::NumericVector prob = Rcpp::NumericVector::create());
arma::uvec sample_cpp(const arma::uvec &x,
                     const int size,
                     const bool replace,
                     Rcpp::NumericVector prob = Rcpp::NumericVector::create());
#endif
