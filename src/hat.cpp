#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List hatcpp(const Rcpp::List& obj){

  const arma::vec& lambda = obj["lambda"];
  const arma::vec& Ed = obj["Ed"];
  const arma::mat& EU = obj["EU"];
  arma::mat carp = obj["carp"];

  int l = lambda.n_elem;
  const arma::mat& rhs = EU.t();
  carp = carp.each_col() % Ed;

  arma::mat ortaveson;
  Rcpp::List hat(l);
  for(int i = 0; i < l; i++){
    ortaveson = rhs.each_col() % carp.col(i);
    hat[i] = EU * ortaveson;
  }

  return hat;
}
