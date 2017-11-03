#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
//' Vague picture
//' 
//' @param pic is the input matrix needs to vague
//' @return the matrix which has already vagued
//' @examples
//' require(recoverpic)
//' vimg = vague(img)
// [[Rcpp::export]]
arma::mat vague(arma::mat pic) {
  int ncol   = pic.n_cols;    // dimension p
  int nrow   = pic.n_rows;    // n
  int noinum = ncol*nrow*0.4;
  int i;
  ivec rran   = randi(noinum, distr_param(0,nrow-1));
  ivec cran   = randi(noinum, distr_param(0,ncol-1));
  for(i=0;i<noinum;i++)
  {
    pic(rran(i),cran(i)) = 100;
  }
  return pic;
}