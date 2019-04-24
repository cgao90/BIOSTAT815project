#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;


arma::mat genBootIndices(int B, int n);       


//' An Rcpp function that calculates the bootstrap distribution of the principal components (PCs) of a low dimensional matrix.
//'
//' @param  n    A (n x n) matrix of sample scores (i.e. the product of singular values and transpose(right singular vectors) 
//' @param  DUt  transpose of UD (matrix of scores, were rows denote individuals, and columns denote measurements in the PC space.)
//' @param  K    the number of PCs to be estimated.
//' @param  B    the number of bootstrap samples   
//' @param  centerSamples  whether each bootstrap sample should be centered before calculating the SVD, 1=Yes, 0=NO(default).
//' @return A list containing left, right singular vectors and singular values of SVD(sample scores)
//' @export
// [[Rcpp::export]]

List bootSVD_LD(int n, arma::mat& DUt, int K, int B, int centerSamples = 0) {
  List dbs(B), Abs(B), Ubs(B); 
  arma::mat bInds = genBootIndices(B, n);
  arma::mat DUtP(n, n); 

  for (int i = 0; i < B; i++) {   // loop through each bootstrap sample
    for (int j = 0; j < n; j++)   // obtain bootstrapped sample from sample scores
      DUtP.col(j) = DUt.col(bInds(i, j)); 

      if (centerSamples == 1) {  // if choose to center the data
        for (int k = 0; k < n; k++) {
          double mean = 0;
          for (int q = 0; q < n; q++)
            mean += DUtP.row(k)[q];
          mean /= n;
          for (int q = 0; q < n; q++)
            DUtP.row(k)[q] -= mean;
        }
      }

      // SVD on input data Y
      arma::mat V1;
      arma::vec d1;
      arma::mat U1;
      arma::svd_econ(V1, d1, U1, DUtP); // economical SVD 

      int nd = d1.size(); // size of singular values
      arma::mat D1(nd, 1);
      D1.col(0) = d1; // vector  -> n x 1 matrix
      dbs(i) = (D1); //

      vector<int> signSwitcher;
      for (int j = 0; j < n; j++)
        V1(j, j) >= 0 ? signSwitcher.push_back(1) : signSwitcher.push_back(-1);
        
      arma::mat tempAbs(n, K), tempUbs(n, K);
      for (int p = 0; p < K; p++) {
        tempAbs.col(p) = signSwitcher[p] * V1.col(p);
        tempUbs.col(p) = signSwitcher[p] * U1.col(p);
      } 
      
      Abs(i) = tempAbs;
      Ubs(i) = tempUbs;
  }

  // save the low dimension SVD results in list
  List bootSVD_LD_output = List::create(Named("dbs") = dbs,
                                        Named("Abs") = Abs,
                                        Named("Ubs") = Ubs);
  return bootSVD_LD_output;   
}  
  



