#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

void quickSort(int q, arma::mat& A, int p, int r);

//' An Rcpp function that calculate bootstrap percentile-based confidence intervals for the PCs of sample socores.
//' @param  VsByK    A list of Vs indexed by K
//' @param  B        number of bootstrap samples
//' @param  K        first K PCs of interest
//' @param  p        number of features
//' @param  alpha    (100*alpha/2, 100(1-alpha/2)) CI
//' @return A list of bootstrap percentile-based confidence intervals for the PCs.
//' @export
// [[Rcpp::export]]

List getLDpercentiles(List& AsByK, int B, int K, int n, double alpha) {  
  List LDpercentiles(K);
  vector<double> percentiles;
  percentiles.push_back(alpha / 2);
  percentiles.push_back(1.0 - alpha / 2);

  for (int k = 0; k < K; k++) {
    arma::mat Ask = AsByK(k); // B x p boostrapped Vs for kth PC 
    Ask = Ask.t();            // p x B
    arma::mat Ak_temp(n, 2);  // to store lb and ub of quantile intervals for kth PC
        //temp = MatrixXf::Zero(p, 2);
    for (int i = 0; i < n; i++) {
      quickSort(i, Ask, 0, B - 1); // sort each row
      double a = percentiles[0] * (B - 1) + 1;
      double b = a - int(a);
      double q1 = double(Ask(i, int(a) - 1)) + double(b * (Ask(i, int(a)) - Ask(i, int(a) - 1)));
      a = percentiles[1] * (B - 1) + 1;
      b = a - int(a);
      double q2 = Ask(i, int(a) - 1) + b * (Ask(i, int(a)) - Ask(i, int(a) - 1));
      Ak_temp.row(i)[0] = q1;
      Ak_temp.row(i)[1] = q2;
    }
    LDpercentiles(k) = Ak_temp;
  } 

  return LDpercentiles;                        
}


