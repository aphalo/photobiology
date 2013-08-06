#include <Rcpp.h>
using namespace Rcpp;

//' Gives irradiance from spectral irradiance.
//'
//' This function gives the result of integrating spectral irradiance over
//' wavelengths. Coded in C++.
//' 
//' @usage integrate_irradianceC(w_length, s_irrad)
//' 
//' @param w_length numeric array of wavelength (nm)
//' @param s_irrad numeric array of spectral irradiances
//' 
//' @return a single numeric value with no change in scale factor: e.g. [W m-2 nm-1] -> [W m-2]
//' @keywords manip misc
//' @name integrate_irradianceC
//' @alias integrate_irradianceC
//' @export
//' @useDynLib photobiology
//' @examples
//' data(sun.data)
//' with(sun.data, integrate_irradianceC(w.length, s.e.irrad))
// [[Rcpp::export]]
double integrate_irradianceC(NumericVector w_length, NumericVector s_irrad) {
    double irradiance = 0.0;
    int n = w_length.size();
    for (int i = 0; i < n-1; i++){
      irradiance =+ ((s_irrad[i+1] + s_irrad[i]) / 2.0) * (w_length[i+1] - w_length[i]);
    }
    return irradiance;
}
