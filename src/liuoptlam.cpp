#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


double bilin_form(const arma::vec &x, const arma::mat &A, const arma::vec &y){
  return sum(x%(A*y));
}

arma::mat inv_A_xxtcpp(const arma::mat &XtX, const arma::vec &x){
  arma::mat XtXinv = arma::inv_sympd(XtX);
  arma::vec XtXinvx = XtXinv * x;
  double karf = sum(x % (XtXinv * x));
  return XtXinv + (XtXinvx * XtXinvx.t()) / (1.0 - karf);
}

// [[Rcpp::export]]
arma::vec liuoptlamcpp(const Rcpp::List& obj){

  arma::vec optlams(5);

  const arma::vec& lambda = obj["lambda"];
  const arma::vec& Ed = obj["Ed"];
  const arma::mat& EU = obj["EU"];
  const arma::mat& EV = obj["EV"];
  const arma::mat& carp = obj["carp"];
  const arma::vec& coefls = obj["coefls"];
  const arma::rowvec& SSRES = obj["SSRES"];
  const arma::vec& residls = obj["residls"];
  const arma::mat& X = obj["X"];
  const arma::vec& y = obj["y"];

  const arma::vec& Ed2 = arma::square(Ed);
  const arma::mat& EUt = EU.t();

  int n = X.n_rows;
  int p = X.n_cols;

  double sigma2 = dot(residls, residls)/(1.0*n - 1.0*p);

  const arma::vec& alphals = EV.t()*coefls;
  const arma::vec& alphals2 = arma::square(alphals);

  // MAKALELERDEN
  const arma::vec& Ed21 = Ed2 + 1.0;

  // //Liu, 1993 lammm
  const arma::vec& numlmm = 1.0 / (Ed2 % Ed21);
  const arma::vec& dnumlmm = (alphals2) / (arma::square(Ed21));
  optlams(0) = 1.0 - sigma2 * (sum(numlmm) / sum(dnumlmm));

  // //Liu, 1993 lamcl
  const arma::vec& numlcl = 1.0 / Ed21;
  const arma::vec& dnumlcl = (Ed2 % alphals2) / (arma::square(Ed21));
  optlams(1) = 1.0 - sigma2 * (sum(numlcl) / sum(dnumlcl));

  // Liu, 1993 lamopt
  optlams(2) = sum((alphals2 - sigma2)/arma::square(Ed21))/
    sum((sigma2 + Ed2 % alphals2)/(Ed2 % arma::square(Ed21)));

  // Improved lam OZKALE
  if(n > p){
    arma::vec ehat = arma::vec(n);
    arma::vec etild = arma::vec(n);
    const arma::mat& Xt = X.t();
    const arma::mat& XtX = X.t() * X;
    arma::mat Xty_Xyt = ((X.each_col() % y).t());
    Xty_Xyt = Xt * y - Xty_Xyt.each_col();
    const arma::vec& diagH = arma::sum(arma::square(EU), 1); // = diag(tcrossprod(EU))  BURADA  H = (X%*%solve(Xtx))%*%t(X)
    const arma::vec& diagG = arma::sum(arma::square(((EUt.each_col() % arma::sqrt(Ed2/Ed21)).t())), 1);
    for(int i = 0; i < n; i++){
      ehat(i) = y(i)-bilin_form(Xt.col(i), inv_A_xxtcpp(XtX, Xt.col(i)), Xty_Xyt.col(i));
      etild(i) = y(i)-bilin_form(Xt.col(i), inv_A_xxtcpp(XtX + arma::diagmat(arma::vec(arma::ones(p))), Xt.col(i)), Xty_Xyt.col(i));
    }
    //
    const arma::vec& F1 = (etild/(1.0-diagG)) % (etild/(1.0-diagH) - ehat/(1.0-diagH));
    const arma::vec& F2 = arma::square(etild/(1.0-diagG) - (ehat/(1.0-diagH)));
    optlams(3) = sum(F1)/sum(F2);
  } else {
    optlams(3) = arma::datum::nan;
  }

  // GCV
  const arma::mat& Ecarpmat = carp.each_col()%Ed;
  const arma::rowvec& GCVliu = SSRES/arma::square(n - 1.0 - arma::sum(Ecarpmat));
  optlams(4) = lambda(GCVliu.index_min());

  return  optlams;
}
