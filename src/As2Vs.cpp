#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

//' An Rcpp function that calculates the bootstapped Vs
//' @param  AsByK    number of bootstrap samples indexed by K
//' @param  V        left singular vectors of original sample data
//' @return A list of bootstrapped Vs indexed by Bootstrap samples
//' @export
// [[Rcpp::export]]

List As2Vs(List& AsByB, arma::mat V) {
  List VsByB(AsByB.size());
  for (int i = 0; i < AsByB.size(); i++) {
    arma::mat A_b = AsByB[i];
    VsByB(i) = V * A_b;
  }
  return VsByB;
}
