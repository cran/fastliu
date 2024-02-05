#ifndef predictliureg_H
#define predictliureg_H

#include <RcppArmadillo.h>

arma::mat predict_liureg(const Rcpp::List& obj, const arma::mat& newdata);

#endif
