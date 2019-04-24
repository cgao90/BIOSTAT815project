#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;



arma::mat genBootIndices(int B, int n);

List bootSVD_LD(int n, arma::mat& DUt, int K, int B, int centerSamples = 0);

List reindexMatricesByK(List& AsByB, int B, int K);

List getMomentsAndMomentCI(List& AsByK, arma::mat& V, int B, int K, int n);

void quickSort(int q, arma::mat& A, int p, int r);

List As2Vs(List& AsByB, arma::mat V);

List getHDpercentiles(List& VsByK, int B, int K, int p, double alpha);

List getLDpercentiles(List& AsByK, int B, int K, int n, double alpha);

//' An Rcpp function that carries out bootstrap PCA and calculats moment-based CIs and percentile CIs
//'
//' @param  Y    A (p x n) matrix of observed data points. 
//' @param  K   num of leading PCs of interest
//' @param  B   number of bootstrap samples
//' @param  alpha  parameter for calculating percentile CIs
//' @return A list of boostrap results
//' @export
// [[Rcpp::export]]

List bootPCA(arma::mat& Y,            // input data (p measurements of n subjects)
             int K,                   // num of leading PCs of interest
             int B,                   // number of bootstrap samples
             double alpha,
             int centerSamples = 0) {  // num of bootstrap samples 

  int p = Y.n_rows;
  int n = Y.n_cols;

  // SVD on input data Y
  arma::mat V; // left singular vectors
  arma::vec D; // singular values
  arma::mat U; // right singular vectors
  arma::svd_econ(V, D, U, Y);

  cout << "SVD on High Dimensional Data Completed" << endl;

  // UD_transpose
  arma::mat DUt = arma::diagmat(D) * U.t();
  //cout << "hello1" << endl;

  //int n = DUt.n_cols;

  // bootSVD_LD(low_dimension)
  List bootSVD_LD_output = bootSVD_LD(n, DUt, K, B, centerSamples);
  cout << "Bootstrap SVD on Low Dimensional Data Completed" << endl;


  // full_LD_PC (Abs)
  List AsByB = bootSVD_LD_output["Abs"]; 
  
  //As vector indexed by K (PC index), each item has dim=BxK
  List AsByK = reindexMatricesByK(AsByB, B, K); 

  //cout << "hello2" << endl;

  // LD moments and moment CI
  int diag1 = min(Y.n_rows, Y.n_cols);
  arma::vec diag_vec = arma::ones(diag1);
  arma::mat diag_mat = diagmat(diag_vec); //generate an identity matrix
  List momentsAndMomentCI_LD = getMomentsAndMomentCI(AsByK, diag_mat, B, K, n);
  cout << "Low Dimensional moments and moment CI Calculated" << endl;


  // HD moments and moment CI
  List momentsAndMomentCI_HD = getMomentsAndMomentCI(AsByK, V, B, K, n);

  cout << "High Dimensional moments and moment CI Calculated" << endl;

  //full_HD_PC_dist
  List VsByB = As2Vs(AsByB, V);   //Vs is full_HD_PC_dist
  //VsByK
  List VsByK= reindexMatricesByK(VsByB, B, K);

  // LD percentile CI
  List percentilesCI_LD = getLDpercentiles(AsByK, B, K, n, alpha);
  cout << "Low Dimensional percentile CI Calculated" << endl;

  
  // HD percentile CI
  List percentilesCI_HD = getHDpercentiles(VsByK, B, K, p, alpha);
  cout << "High Dimensional percentile CI Calculated" << endl;

  


  List bootPCA_output = List::create(Named("V") = V,
                                     Named("bootSVD_LD_output") = bootSVD_LD_output,
                                     Named("momentsAndMomentCI_LD") = momentsAndMomentCI_LD,
                                     Named("momentsAndMomentCI_HD") = momentsAndMomentCI_HD,
                                     Named("full_LD_PC_dist") = AsByK,
                                     Named("full_HD_PC_dist") = VsByK,
                                     Named("percentilesCI_LD") = percentilesCI_LD,
                                     Named("percentilesCI_HD") = percentilesCI_HD);



  return bootPCA_output;

}



