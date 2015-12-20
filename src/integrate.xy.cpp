#include <Rcpp.h>
using namespace Rcpp;

using namespace std;

//' Gives irradiance from spectral irradiance.
//'
//' This function gives the result of integrating spectral irradiance over
//' wavelengths. Coded in C++.
//'
//' @param x numeric array
//' @param y numeric array
//'
//' @return a single numeric value with no change in scale factor: e.g. [W m-2
//'   nm-1] -> [W m-2]
//' @keywords manip misc
//' @export
//' @useDynLib photobiology
//' @examples
//' with(sun.data, integrate_xy(w.length, s.e.irrad))
// [[Rcpp::export]]
double integrate_xy(NumericVector x, NumericVector y) {
    int n = x.size();
    double irradiance = 0.0;
    for (int i = 0; i < n-1; i++){
      irradiance += ((y[i+1] + y[i]) / 2.0) * (x[i+1] - x[i]);
    }
    return irradiance;
}
