#ifndef liureg_H
#define liureg_H

#include <RcppArmadillo.h>

Rcpp::List liuregcpp(const arma::mat &Xmat,
                     const arma::vec &yvec,
                     const arma::vec &lambda,
                     std::string scale);

#endif
