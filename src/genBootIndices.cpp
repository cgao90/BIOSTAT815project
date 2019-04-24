#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

//' An Rcpp function that generate bootstrapped indices
//' @param  B    number of bootstrap samples
//' @param  n    number of original sample
//' @return A (B x n) matrix  containing the bootstrapped indices
//' @export
// [[Rcpp::export]]

arma::mat genBootIndices(int B,
                         int n) {       
  arma::mat bInds(B, n);
  for (int i = 0; i < B; i++) {
    for (int j = 0; j < n; j++) {
      bInds(i, j) = rand() % n;
    }
  }
  return bInds;                       
}


