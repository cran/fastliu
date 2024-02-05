// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "coefliureg.h"

// [[Rcpp::export]]
arma::mat predict_liureg(const Rcpp::List& obj, const arma::mat& newdata) {
  int nn = newdata.n_rows;
  const arma::mat& beta = coef_liureg(obj);
  const arma::vec& bir = arma::ones(nn);
  arma::mat newdatafit = arma::join_rows(bir, newdata) * beta;
  return newdatafit;
}
