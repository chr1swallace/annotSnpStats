#include<Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
RawMatrix asw(RawMatrix x, NumericVector w) {
  int nrow=x.nrow(),
    nw=w.size();
  RawMatrix nx = clone(x);
  // cout << "hello";
  for(int jw=0; jw<nw; jw++) {
    int j=w[jw]-1;  // 1-based to 0-based
    for(int i=0; i<nrow; i++)
      if(nx(i,j)==1 || nx(i,j)==3)
	nx(i,j) = 4 - nx(i,j);
  }
  return nx;
}

