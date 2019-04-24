#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

//' An Rcpp function to order values in each row of the input matrix
//'
//' @param  q    number of dimension in features
//' @param  A    A matrix to be ordered
//' @param  p 
//' @param  r 

//' @return A matrix whose values in each row has ascending order
//' @export
// [[Rcpp::export]]

void quickSort(int q, arma::mat& A, int p, int r) {
  if (p < r) { // immediately terminate if subarray size is 1
    double piv = A(q, r); // take a pivot value
    int i = p - 1; // p-i-1 is the # elements < piv among A[p..j]
    double tmp;
    for (int j = p; j < r; j++) {
      if (A(q, j) < piv) { // if smaller value is found, increase q (=i+1)
        i++;
        tmp = A(q, i); A(q, i) = A(q, j); A(q, j) = tmp; // swap A[i] and A[j]
      }
    }
    A(q, r) = A(q, i + 1); A(q, i + 1) = piv; // swap A[i+1] and A[r]
    quickSort(q, A, p, i);
    quickSort(q, A, i + 2, r);
  }
}

