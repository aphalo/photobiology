#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;

  for(int i = 0; i < n; ++i) {
    total += x[i] / n;
  }
  return total;
}

/*** R 
  library(microbenchmark)
  x <- runif(1e5)
  microbenchmark(
    mean(x),
    meanC(x))
*/