#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List liuregcpp(const arma::mat &Xmat,
                     const arma::vec &yvec,
                     const arma::vec &lambda,
                     std::string scale) {

  int n = Xmat.n_rows;
  int p = Xmat.n_cols;

  const arma::rowvec& Xm = arma::mean(Xmat);
  double ym = mean(yvec);

  arma::mat X = Xmat.each_row()-Xm;
  arma::vec y = yvec - ym;

  arma::rowvec Xs;
  if(scale == "ulength"){
    Xs = arma::sqrt(arma::sum(arma::square(X)));
  } else if(scale == "unormal"){
    Xs = arma::sqrt(arma::sum(arma::square(X))/double(n));
  } else{
    Xs = arma::ones<arma::rowvec>(p);
  }
  X.each_row() /= Xs;

  int np = std::min(n,p);
  arma::mat EU;
  arma::mat EV;
  arma::vec Ed;
  arma::svd_econ(EU, Ed, EV, X, "both", "dc");
  EU = EU.head_cols(np);
  EV = EV.head_cols(np);

  const arma::vec& Ed2 = arma::square(Ed);
  const arma::vec& Ed12 = (Ed2 + 1.0)%Ed;

  const arma::vec& rhs = EU.t() * y;

  const arma::vec& coefls = EV * (rhs/Ed);
  const arma::vec& fitls = X * coefls;
  const arma::vec& residls = y - fitls;

  int l = lambda.n_elem;
  int dx = Ed.n_elem;
  arma::mat carp = arma::repmat(lambda,1,dx).t();
  for(int i = 0; i < l; i++){
    arma::vec carpptr(carp.colptr(i),dx,false,true);
    carp.col(i) = (carpptr + Ed2)/Ed12;
  }
  const arma::mat& ortaveson = rhs%carp.each_col();
  // carp.each_col() += Ed2;
  // carp.each_col() %= rhs/(Ed12);
  const arma::mat& coefliu = EV * ortaveson;

  const arma::mat& fitliu = X * coefliu;
  const arma::mat& residliu = y - fitliu.each_col();
  //
  double SST = arma::dot(y,y);
  const arma::rowvec& SSRES = arma::sum(arma::square(residliu));

  return Rcpp::List::create(Rcpp::Named("coefliu") = coefliu,
                            Rcpp::Named("fitliu") = fitliu,
                            Rcpp::Named("residliu") = residliu,
                            Rcpp::Named("Ed") = Ed,
                            Rcpp::Named("EU") = EU,
                            Rcpp::Named("EV") = EV,
                            Rcpp::Named("carp") = carp,
                            Rcpp::Named("SST") = SST,
                            Rcpp::Named("SSRES") = SSRES,
                            Rcpp::Named("coefls") = coefls,
                            Rcpp::Named("fitls") = fitls,
                            Rcpp::Named("residls") = residls,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("X") = X,
                            Rcpp::Named("y") = y,
                            Rcpp::Named("Xscale") = Xs,
                            Rcpp::Named("Xm") = Xm,
                            Rcpp::Named("ym") = ym,
                            Rcpp::Named("scale") = scale);
}
