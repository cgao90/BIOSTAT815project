#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

//' An Rcpp function that calculate bootstrap moments and moment-based confidence intervals for the PCs.
//'
//' @param  AsByK  number of bootstrap samples indexed by K
//' @param  V      left singular vectors of original sample data
//' @param  B      number of bootstrap samples
//' @param  K      first K PCs of interest
//' @param  n      number of subjects 
//' @return A list of moments and moment-based CIs for the PCs
//' @export
// [[Rcpp::export]]

List getMomentsAndMomentCI(List& AsByK, arma::mat& V, int B, int K, int n) {    
  int p = V.n_rows; //number of features
  arma::vec ones = arma::ones(B); // column vector of ones
  List EPCs(K), varVs(K), sdVs(K), momentCI(K);  // expectation/variance/standard deviation/95% CI of PCs
  
  for (int i = 0; i < K; i++) {
    // EAs: 1*n matrixs
    arma::mat AsByK_temp = AsByK[i];
    arma::mat EAs = arma::mean(AsByK_temp); // column means 1 x n
    // EPCs: p*n * n*1 = p*1
    arma::mat EPCs_temp = V * EAs.t();
    EPCs(i) = EPCs_temp;
    
    // Ones: n*1
    // varAs: p*n * n*B * B*n
    // varVs
    arma::mat varAs = V * ((AsByK_temp - ones * EAs).t() * (AsByK_temp - ones * EAs)) / (n - 1);
    
    // temp: p*n * p*n = p*n
    arma::mat varVs_temp = varAs % V;
    varVs_temp = sum(varVs_temp, 1);
    varVs(i) = varVs_temp;

    // sdVs
    arma::mat sdV_temp = arma::sqrt(varVs_temp);
    sdVs(i) = sdV_temp;
    // CI
    arma::mat momentCI_temp(p, 2);
    momentCI_temp.col(0) = EPCs_temp - 1.96 * sdV_temp;
    momentCI_temp.col(1) = EPCs_temp + 1.96 * sdV_temp;
    momentCI(i) = momentCI_temp;
  }
  // store moments and CI in vector
  List MomentsAndMomentCI_output = List::create(Named("EPCs") = EPCs,
                                                Named("varVs") = varVs,
                                                Named("sdVs") = sdVs,
                                                Named("momentCI") = momentCI);

  return MomentsAndMomentCI_output;                        
}


