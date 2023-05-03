#include<algorithm>
#include<RcppArmadillo.h>
#include<omp.h>
using namespace Rcpp;
using namespace std;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List COMPSSE(arma::mat Xo, arma::mat Yo, arma::vec Z, double r1, double r2, int N) {

	 arma::mat BETAo,SSE,YB,M1,M2,M3;
   int I1,I2,COUNT;

  for(int L=0; L < N ; L++){
        if(r1 < Z[L] && Z[L]<=r2){
      	  I1=1;
        }else{I1=0;}
      if(L==0){COUNT=I1;}else{COUNT=COUNT+I1;}
      YB=Xo(span(L),span::all).t();
         if(L==0){
          M1=(YB*YB.t())*I1;
          M2=(YB*Yo(span(L),span::all))*I1;
         }else{
          M1=M1+((YB*YB.t())*I1);	
          M2=M2+((YB*Yo(span(L),span::all))*I1);
         }
     }
      BETAo=pinv(M1,0.01)*M2;

    for(int M=0; M < N ; M++){
        if(r1 < Z[M] && Z[M]<=r2){
          I2=1;
        }else{I2=0;}
      YB=Xo(span(M),span::all).t();
         if(M==0){
           M3=pow(((Yo(span(M),span::all)-(YB.t()*BETAo))*I2),2);
         }else{
           M3=M3+pow(((Yo(span(M),span::all)-(YB.t()*BETAo))*I2),2);
         }
     }
      SSE=M3;

return Rcpp::List::create(
    Rcpp::Named("BETA")=BETAo,
    Rcpp::Named("SSE")=SSE,
    Rcpp::Named("COUNT")=COUNT
    );
}