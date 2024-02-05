#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat liustatscpp(const Rcpp::List& obj){

  const arma::vec& lambda = obj["lambda"];
  const arma::mat& coefliu = obj["coefliu"];
  const arma::mat& fitliu = obj["fitliu"];
  const arma::mat& residliu = obj["residliu"];
  const arma::vec& Ed = obj["Ed"];
  // const arma::mat& EU = obj["EU"];
  const arma::mat& EV = obj["EV"];
  const arma::mat& carp = obj["carp"];
  double SST = obj["SST"];
  const arma::rowvec& SSRES = obj["SSRES"];
  const arma::vec& coefls = obj["coefls"];

  const arma::rowvec& SSR = arma::sum(arma::square(fitliu));

  const arma::vec& Ed2 = arma::square(Ed);
  const arma::vec& EdEd21 = Ed % (Ed2 + 1.0);

  const arma::mat& EVt = EV.t();

  int n = fitliu.n_rows;
  int p = coefliu.n_rows;
  int l = lambda.n_elem;

  const arma::mat& Ecarpmat = carp.each_col()%Ed;

  const arma::rowvec& edfliu = n - arma::sum(2.0*Ecarpmat-arma::square(Ecarpmat));

  const arma::rowvec& sigma2liu = SSRES/edfliu;

  arma::vec varliu(l);
  arma::mat covliu(p, p*l);
  for(int j = 0; j < l; j++){
    arma::vec covdd = arma::square(((Ed2 + lambda(j))/EdEd21));
    covliu.cols(p*j, p*(j+1)-1) = sigma2liu(j)*(EV * (covdd % EVt.each_col()));
    varliu(j) = sum(arma::diagvec(covliu.cols(p*j, p*(j+1)-1)));
  }

  const arma::vec& abeta = EVt*coefls;
  const arma::vec& bias2liu = arma::square((lambda - 1.0)) *
    sum((arma::square(abeta))/arma::square(Ed2 + 1.0));

  const arma::vec& mseliu = varliu + bias2liu;

  arma::vec Fvliu = arma::vec(l);
  Fvliu.fill(arma::datum::nan);
  if(n > p){
    for(int j = 0; j < l; j++){
      Fvliu(j) = sum(coefliu.col(j)%(arma::solve(covliu.cols(p*j, p*(j+1)-1), coefliu.col(j))))/(p*1.0);
    }
  }

  const arma::rowvec& GCVliu = SSRES/arma::square(n - 1.0 - arma::sum(Ecarpmat));

  const arma::rowvec& R2liu = 1.0-SSRES/SST;

  const arma::rowvec& adjR2liu = 1.0-(n-1.0)/(n-p-1.0)*(1.0-R2liu);

  arma::mat stats(l,9);
  stats.col(0) = edfliu.t();
  stats.col(1) = sigma2liu.t();
  stats.col(2) = varliu;
  stats.col(3) = bias2liu;
  stats.col(4) = mseliu;
  stats.col(5) = Fvliu;
  stats.col(6) = GCVliu.t();
  stats.col(7) = R2liu.t();
  stats.col(8) = adjR2liu.t();

  return stats;

  // return Rcpp::List::create(Rcpp::Named("EDF") = edfliu.t(),
  //                           Rcpp::Named("sigma2") = sigma2liu.t(),
  //                           Rcpp::Named("VAR") = varliu,
  //                           Rcpp::Named("BIAS2") = bias2liu,
  //                           Rcpp::Named("MSE") = mseliu,
  //                           Rcpp::Named("FVal") = Fvliu,
  //                           Rcpp::Named("GCV") = GCVliu.t(),
  //                           Rcpp::Named("R2") = R2liu.t(),
  //                           Rcpp::Named("AdjR2") = adjR2liu.t());
}
