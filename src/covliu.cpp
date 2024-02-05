#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List covliucpp(const Rcpp::List& obj){

  const arma::vec& lambda = obj["lambda"];
  const arma::vec& Ed = obj["Ed"];
  const arma::mat& EU = obj["EU"];
  const arma::mat& EV = obj["EV"];
  const arma::mat& carp = obj["carp"];
  const arma::rowvec& SSRES = obj["SSRES"];

  const arma::vec& Ed2 = arma::square(Ed);
  const arma::vec& EdEd21 = Ed % (Ed2 + 1.0);

  const arma::mat& EVt = EV.t();

  int n = EU.n_rows;
  int l = lambda.n_elem;


  const arma::mat& Ecarpmat = carp.each_col() % Ed;
  const arma::rowvec& edfliu = n - arma::sum(2.0*Ecarpmat-arma::square(Ecarpmat));

  const arma::rowvec& sigma2liu = SSRES/edfliu;

  Rcpp::List covliu(l);
  for(int j = 0; j < l; j++){
    arma::vec covdd = arma::square(((Ed2 + lambda(j))/EdEd21));
    covliu[j] = sigma2liu(j)*(EV * (covdd % EVt.each_col()));
  }

  return covliu;
}
