#ifndef __gen_Ft__
#define __gen_Ft__

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
// [[Rcpp::export]]
arma::mat gen_Ft(arma::mat F1) {
  int n = F1.n_rows;
  arma::mat temp_I;
  temp_I.eye(n, n);
	arma::mat F1t = arma::kron(F1, temp_I);
	return F1t;
}

#endif // __gen_Ft__
