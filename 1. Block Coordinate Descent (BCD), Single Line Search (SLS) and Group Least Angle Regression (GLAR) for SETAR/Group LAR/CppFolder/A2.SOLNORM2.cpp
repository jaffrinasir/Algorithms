#include<algorithm>
#include<RcppArmadillo.h>
#include<omp.h>
using namespace Rcpp;
using namespace std;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List SOLNORM2(arma::mat DATAC, arma::mat resi, arma::mat C0, arma::vec jT, int core) {
  int sizeA = jT.size();
  int sizeB = DATAC.n_rows;
  arma::vec aj(sizeA);
  arma::vec bj(sizeA);
  arma::vec cj(sizeA);
  omp_set_num_threads(core);

  # pragma omp parallel
  {
  arma::mat Ba,Bb;
  # pragma omp for
  for(int L=0; L < sizeA ; L++){
     for(int M=(jT[L]-1); M < sizeB ; M++){  
          if(M == (jT[L]-1)){
            Ba=DATAC(span(M),span::all).t()*C0(span(M),span::all);
            Bb=DATAC(span(M),span::all).t()*resi(span(M),span::all);
          }else{
            Ba=Ba+DATAC(span(M),span::all).t()*C0(span(M),span::all);            
            Bb=Bb+DATAC(span(M),span::all).t()*resi(span(M),span::all);
          }
     }
         bj[L]=sum(vectorise(Ba.t()*Bb));
         aj[L]=sum(vectorise(pow(Bb,2)));
         cj[L]=sum(vectorise(pow(Ba,2)));
  }
 } 
    return Rcpp::List::create(
    Rcpp::Named("aj")=aj,
    Rcpp::Named("bj")=bj,
    Rcpp::Named("cj")=cj
    );
}