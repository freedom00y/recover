#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
arma::mat pOmega(vec pos, arma::mat X);
arma::mat svth(arma::mat A, double lambda);
double fnorm(arma::mat X);

//' Fill in the missing picture
//' 
//' @param pic is the input matrix which needs to be full of
//' @param lambda is the penalty variable
//' @param ac is the bool variable, 1 means using accelerate method, 0 means no accelerate
//' @return the matrix which has already full of
//' @examples
//' require(recoverpic)
//' nimgac = recover(vimg,0.3)
// [[Rcpp::export]]
Rcpp::List recover(arma::mat &pic,double lambda,bool ac=1)
{
  int ncol = pic.n_cols;    // dimension p
  int nrow = pic.n_rows;    // n
  int null = 0, count=0;
  double tol=0.01;
  double t0=1,t1=1;
  vec pos(ncol*nrow);//record Omega position
  mat Y(nrow, ncol, fill::zeros);
  mat hatY(nrow, ncol, fill::zeros);
  mat X=pic, decopo;
  for(int i=0;i<nrow*ncol;i++)
    if(!(pic[i]==NA))
    {
      pos(null)=i;
      null++;
    }
    
    do
    {
      Y=hatY;
      decopo=Y+pOmega(pos.head(null),X-Y);
      hatY=svth(decopo,lambda/2);
      if(ac)
      {
        t0=t1;
        t1=0.5*(1+sqrt(1+4*t0*t0));
        hatY=hatY+(t0-1)/t1*(hatY-Y);
      }
      count++;
    } while (fnorm(Y-hatY)/fnorm(Y)>tol);
    
    for(int i=0;i<nrow*ncol;i++)
    {
      if(hatY(i)<0)
        hatY(i)=0;
      else if(hatY(i)>1)
        hatY(i)=1;
    }
    
    //return hatY;
    return List::create(Named("newimg") = hatY,
                        Named("repeat time")= count
    );
}

arma::mat pOmega(vec pos, arma::mat X)
{
  int ncol = X.n_cols;    // dimension p
  int nrow = X.n_rows;    // n
  int n = pos.n_elem;
  mat Y(nrow, ncol, fill::zeros);
  for(int i=0;i<n;i++)
  {
    Y[pos(i)] = X[pos(i)];
  }
  return Y;
}

arma::mat svth(arma::mat A, double lambda)
{
  int n=A.n_rows;
  int p=A.n_cols;
  mat U,V,B;//B for return
  mat D(n,p,fill::zeros);
  vec d;
  svd(U,d,V,A);
  for(int i=0;i<d.n_elem;i++)
  {
    vec vv(2);
    vv(0)=d(i)-lambda;
    vv(1)=0;
    D(i,i)=max(vv);
  }
  B=U*D*V.t();
  return(B);
}

double fnorm(arma::mat X)
{
  return sqrt(trace(X.t()*X));
}