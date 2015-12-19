#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//' Insert interpolated spectral irradiance values.
//'
//' Inserting wavelengths values immediately before and after a discontinuity in
//' the SWF, greatly reduces the errors caused by interpolating the weighted
//' irradiance during integration of the effective spectral irradiance. This is
//' specially true when data has a large wavelength step size.
//'
//' @param x numeric array of wavelength (nm)
//' @param y numeric array of spectral irradiance values
//' @param hinges a numeric array giving the wavelengths at which the s.irrad
//'   should be inserted by interpolation, no interpolation is indicated by an
//'   empty array (numeric(0))
//'
//' @return a numeric vector
//' @name put_hinges
//' @keywords manip misc
//' @export
//'
// I have tried to optimize this function as much as possible.
// It assumes that hinges are sorted in increasing order.
// It uses binary search to find the insertion point, and also does this search in the x
// vector only in the "tail" starting at the previously inserted hinge.
// The second example tests several boundary conditions: response for hinges outside the range of the spectrum, two succesive
// insertions between the a single pair of values in the spectral data, and for insertion of
// hinges at walengths already present in the data. A test for the boundary condition
// when an empty numeric vector is passed to hinges fails.
// with(sun.data, insert_hinges(w.length, s.e.irrad, c())) in Rcpp
//
// [[Rcpp::export]]
NumericVector put_hinges(NumericVector x, NumericVector y, NumericVector hinges) {

    vector<double> x_in = Rcpp::as< vector<double> >(x);
    vector<double> y_in = Rcpp::as< vector<double> >(y);
    vector<double> hinges_in = Rcpp::as< vector<double> >(hinges);

    vector<double>::iterator h_low = lower_bound(hinges_in.begin(), hinges_in.end(), x_in.front());
    vector<double>::iterator h_high = lower_bound(hinges_in.begin(), hinges_in.end(), x_in.back());

    vector<double>::iterator x_in_it = x_in.begin();
    vector<double>::iterator y_in_it = y_in.begin();

    vector<double>::size_type max_len = x_in.size() + hinges.size();
    vector<double> x_out (max_len); // x_out.reserve(max_len);
    vector<double> y_out (max_len); //y_out.reserve(max_len);

    vector<double>::iterator x_out_it = x_out.begin();
    vector<double>::iterator y_out_it = y_out.begin();

    vector<double>::iterator wl_next_chunk_begin_it;
    vector<double>::iterator si_next_chunk_begin_it = y_in_it;

    for (vector<double>::iterator it=h_low; it<h_high; it++) {
      wl_next_chunk_begin_it = lower_bound(x_in_it-1, x_in.end(), *it);
      vector<double>::size_type num_elements = distance(x_in_it, wl_next_chunk_begin_it);
      si_next_chunk_begin_it = si_next_chunk_begin_it + num_elements;
      x_out_it = copy(x_in_it, wl_next_chunk_begin_it, x_out_it);
      y_out_it = copy(y_in_it, si_next_chunk_begin_it, y_out_it);
      if (*it < *wl_next_chunk_begin_it) {
        *x_out_it = *it;
         double interpol_irrad = *(si_next_chunk_begin_it-1) +
               (*si_next_chunk_begin_it - *(si_next_chunk_begin_it-1)) /
               (*wl_next_chunk_begin_it - *(wl_next_chunk_begin_it-1)) * (*it - *(wl_next_chunk_begin_it-1));
        *y_out_it = interpol_irrad;
        ++x_out_it;
        ++y_out_it;
      }
    x_in_it = wl_next_chunk_begin_it;
    y_in_it = si_next_chunk_begin_it;
    }
    x_out_it = copy(x_in_it, x_in.end(), x_out_it);
    y_out_it = copy(y_in_it, y_in.end(), y_out_it);
    x_out.erase(x_out_it,  x_out.end());
    y_out.erase(y_out_it,  y_out.end());

    return Rcpp::wrap(y_out);
}
