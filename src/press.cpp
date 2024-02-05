// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


// [[Rcpp::export]]
arma::vec pressliucpp(const Rcpp::List& obj){

  const arma::vec& lambda = obj["lambda"];
  const arma::vec& Ed = obj["Ed"];
  const arma::mat& EU = obj["EU"];
  const arma::mat& carp = obj["carp"];
  const arma::mat& residliu = obj["residliu"];
  const arma::vec& residls = obj["residls"];

  int n = EU.n_rows;
  int l = lambda.n_elem;

  const arma::mat& tEU2 = arma::square(EU.t());
  const arma::vec& Ed2 = arma::square(Ed);

  arma::mat diagH(n,l);
  arma::vec carpc;
  for(int i = 0; i < l; i++){
    diagH.col(i) = arma::sum(tEU2.each_col() % (carp.col(i) % Ed)).as_col();
  }

  arma::vec olsH = arma::sum(tEU2).as_col();
  arma::vec ridgeH_1 = arma::sum(tEU2.each_col() % (Ed2/(Ed2+1.0))).as_col();
  arma::mat T1 = residliu.each_col()/(1.0-ridgeH_1);
  arma::vec T2 = residls/((1.0 - ridgeH_1)%(1.0 - olsH));
  arma::mat T3 = ridgeH_1 - diagH.each_col();
  T3.each_col() %= T2;

  arma::mat pressliu_mat = T1 - T3;

  arma::rowvec pressliu = arma::sum(arma::square(pressliu_mat));

  return pressliu.t();
}

