#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat diagHcpp(const Rcpp::List& obj){

  const arma::vec& lambda = obj["lambda"];
  const arma::vec& Ed = obj["Ed"];
  const arma::mat& EU = obj["EU"];
  const arma::mat& carp = obj["carp"];

  int n = EU.n_rows;
  int l = lambda.n_elem;

  const arma::mat& tEU2 = arma::square(EU.t());

  arma::mat diagH(n,l);
  for(int i = 0; i < l; i++){
    diagH.col(i) = arma::sum(tEU2.each_col() % (carp.col(i) % Ed)).as_col();
  }

  return diagH;
}
