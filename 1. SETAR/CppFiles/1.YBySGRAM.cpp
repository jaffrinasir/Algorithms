#include<algorithm>
#include<RcppArmadillo.h>
#include<omp.h>
using namespace Rcpp;
using namespace std;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat YBy(arma::mat X, arma::mat Y, int INDEX, int N) {
	arma::mat AA;
	for(int i1 = (INDEX-1); i1 < N ; i1++){
      if(i1==(INDEX-1)){ AA=(X(span(i1),span::all).t())*Y(span(i1),span::all);
      }else{
      	AA=AA+( (X(span(i1),span::all).t())*Y(span(i1),span::all));
      }
	}
	return(AA);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat SGRAM(arma::mat X, int INDEX, int N) {
	arma::mat AB;
	for(int i2 = (INDEX-1); i2 < N ; i2++){
      if(i2==(INDEX-1)){ AB=(X(span(i2),span::all).t())*X(span(i2),span::all);
      }else{
      	AB=AB+( (X(span(i2),span::all).t())*X(span(i2),span::all));
      }
	}
	return(AB);
}