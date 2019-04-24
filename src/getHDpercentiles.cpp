#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

void quickSort(int q, arma::mat& A, int p, int r);

//' An Rcpp function that calculate bootstrap percentile-based confidence intervals for the PCs of original data.
//' @param  VsByK    A list of Vs indexed by K
//' @param  B        number of bootstrap samples
//' @param  K        first K PCs of interest
//' @param  p        number of features
//' @param  alpha    (100*alpha/2, 100(1-alpha/2)) CI
//' @return A list of bootstrap percentile-based confidence intervals for the PCs.
//' @export
// [[Rcpp::export]]

List getHDpercentiles(List& VsByK, int B, int K, int p, double alpha) {  
  List HDpercentiles(K);
  vector<double> percentiles;
  percentiles.push_back(alpha / 2);
  percentiles.push_back(1.0 - alpha / 2);

  for (int k = 0; k < K; k++) {
    arma::mat Vsk = VsByK(k); // B x p boostrapped Vs for kth PC 
    Vsk = Vsk.t();            // p x B
    arma::mat Vk_temp(p, 2);  // to store lb and ub of quantile intervals for kth PC
        //temp = MatrixXf::Zero(p, 2);
    for (int i = 0; i < p; i++) {
      quickSort(i, Vsk, 0, B - 1); // sort each row
      double a = percentiles[0] * (B - 1) + 1;
      double b = a - int(a);
      double q1 = double(Vsk(i, int(a) - 1)) + double(b * (Vsk(i, int(a)) - Vsk(i, int(a) - 1)));
      a = percentiles[1] * (B - 1) + 1;
      b = a - int(a);
      double q2 = Vsk(i, int(a) - 1) + b * (Vsk(i, int(a)) - Vsk(i, int(a) - 1));
      Vk_temp.row(i)[0] = q1;
      Vk_temp.row(i)[1] = q2;
    }
    HDpercentiles(k) = Vk_temp;
  } 

  return HDpercentiles;                        
}


