#include<algorithm>
#include<RcppArmadillo.h>
#include<omp.h>
using namespace Rcpp;
using namespace std;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List QUADcpp(arma::vec aj, arma::vec bj, arma::vec cj, double dj, int core) {
  int sizeA = aj.size();
  arma::vec aPlus(sizeA);
  arma::vec aMinus(sizeA);
  omp_set_num_threads(core);

  # pragma omp parallel
  {
  # pragma omp for
  for(int L=0; L < sizeA ; L++){
         if((cj[L]-dj)==0){
          aPlus[L]=(aj[L]-dj)/2*(bj[L]-dj);
          aMinus[L]=(aj[L]-dj)/2*(bj[L]-dj);
         }else if( (pow((bj[L]-dj),2)-((aj[L]-dj)*(cj[L]-dj))) < 0 ){
          aPlus[L]=1;
          aMinus[L]=1;
         } else{
           aPlus[L]=((bj[L]-dj)+pow((pow((bj[L]-dj),2)-((aj[L]-dj)*(cj[L]-dj))),0.5) )/(cj[L]-dj);
           aMinus[L]=((bj[L]-dj)-pow((pow((bj[L]-dj),2)-((aj[L]-dj)*(cj[L]-dj))),0.5) )/(cj[L]-dj);
         }
   }
 } 
    return Rcpp::List::create(
    Rcpp::Named("aPlus")=aPlus,
    Rcpp::Named("aMinus")=aMinus
    );
}