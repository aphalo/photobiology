#include <Rcpp.h>
using namespace Rcpp;

//' Gives irradiance from spectral irradiance.
//'
//' This function gives the result of integrating spectral irradiance over
//' wavelengths. Coded in C++.
//' 
//' @usage integrateirr(w_length, s_irrad)
//' 
//' @param w_length numeric array of wavelength (nm)
//' @param s_irrad numeric array of spectral irradiances
//' 
//' @return a single numeric value with no change in scale factor: e.g. [W m-2 nm-1] -> [W m-2]
//' @keywords manip misc
//' @name integrateirr
//' @alias integrateirr
//' @export
//' @examples
//' data(sun.data)
//' with(sun.data, integrateirr(w.length, s.e.irrad))
// [[Rcpp::export]]
double integrateirr(NumericVector w_length, NumericVector s_irrad) {
    double irrad = 0.0;
    int last = w_length.size();
    for (int i = 0; i < last; i++){
      irrad =+ ((s_irrad[i] + s_irrad[i+1]) / 2.0 * (w_length[i+1] - w_length[i]));
    }
    return(irrad);
}

