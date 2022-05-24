#include<algorithm>
#include<RcppArmadillo.h>
#include<omp.h>
using namespace Rcpp;
using namespace std;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec CORRFUNCcpp(arma::mat DATAC, arma::mat resi, arma::vec group, int core) {
  int sizeA = group.size();
  int sizeB = DATAC.n_rows;
  arma::vec corr0(sizeA);
  omp_set_num_threads(core);

  # pragma omp parallel
  {
  arma::mat SUMT,SUMT2;
  arma::vec SUMT3;
  # pragma omp for
  for(int L=0; L < sizeA ; L++){
     for(int M=L; M < sizeB ; M++){  
          if(M == L){
            SUMT=DATAC(span(M),span::all).t()*resi(span(M),span::all);
          }else{
            SUMT=SUMT+DATAC(span(M),span::all).t()*resi(span(M),span::all);
          }
            SUMT2=pow(SUMT,2);
            SUMT3=vectorise(SUMT2);
            corr0[L]=pow(sum(SUMT3),0.5);
     }
  }
 } 
  return(corr0); 
}