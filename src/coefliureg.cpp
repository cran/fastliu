// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "liureg.h"

// [[Rcpp::export]]
arma::mat coef_liureg(const Rcpp::List& obj){

  const arma::mat& coefliu = obj["coefliu"];
  const arma::rowvec& Xscale = obj["Xscale"];
  const arma::rowvec& Xm = obj["Xm"];
  double ym = obj["ym"];
  const arma::mat& borj = coefliu.each_col()/Xscale.t();
  const arma::rowvec& b0orj = ym - Xm * borj;
  const arma::mat& betaorj = arma::join_cols(b0orj,borj);
  return betaorj;
}
