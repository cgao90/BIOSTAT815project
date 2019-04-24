#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

//' An Rcpp function that Used for calculation of low dimensional standard errors & percentiles, by re-indexing the \eqn{A^b} by PC index (\eqn{k}) rather than bootstrap index (\eqn{b}).
//'
//' @param  AsByB   a size B list of matrices in the order boostrap samples
//' @param  B       number of boostrap samples
//' @param  K       first K PCs of interest
//' @return A list of reordered matrices (order by PCs)
//' @export
// [[Rcpp::export]]

List reindexMatricesByK(List& AsByB, int B, int K) {     // n = number 
  arma::mat AsByB_eg = AsByB[0];
  int L = AsByB_eg.n_rows; // number of row in each bootstrap sample
  List AsByK(K);

  for (int i = 0; i < K; i++) {
        arma::mat AsByK_temp(B, L); 
        for (int j = 0; j < B; j++) {
          arma::mat AsByB_temp = AsByB[j];
          AsByK_temp.row(j) = AsByB_temp.col(i).t();
        }
        AsByK(i) = AsByK_temp;
  }
  return AsByK;                         // return As reindexed by K
}
