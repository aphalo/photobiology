#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#undef DEBUG
#ifdef DEBUG
#include <iostream>
#endif

//' Insert new wavelengths into the data, interpolating the corresponding spectral irradiance values.
//'
//' Inserting wavelengths values immediately bafore and after a discontinuity in the SWF,
//' gretaly reduces the errors caused by interpolating the weighted irradiance during
//' integration of the effective spectral irradiance. This is specially true when data
//' has a large wavelength step size.
//'
//' @usage insert_hinges(w_length, s_irrad, hinges)
//'
//' @param w_length numeric array of wavelength (nm)
//' @param s_irrad numeric array of spectral irradiance values
//' @param hinges a numeric array giving the wavelengths at which the s.irrad should be inserted by interpolation, no interpolation is indicated by an empty array (numeric(0))
//'
//' @return a data.frame with variables \code{w.length} and \code{s.irrad}
//' @name insert_hinges
//' @keywords manip misc
//' @export
//' @examples
//' with(sun.data, insert_hinges(w.length, s.e.irrad, c(399.99,400.00,699.99,700.00)))
//' with(sun.data, insert_hinges(w.length, s.e.irrad, c(199.99,200.00,399.50,399.99,400.00,699.99,700.00,799.99,1000.00)))
//
// I have tried to optimize this function as much as possible.
// It assumes that hinges are sorted in increasing order.
// It uses binary search to find the insertion point, and also does this search in the w_length
// vector only in the "tail" starting at the previously inserted hinge.
// The second example tests several boundary conditions: response for hinges outside the range of the spectrum, two succesive
// insertions between the a single pair of values in the spectral data, and for insertion of
// hinges at walengths already present in the data. A test for the boundary condition
// when an empty numeric vector is passed to hinges fails.
// with(sun.data, insert_hinges(w.length, s.e.irrad, c())) in Rcpp
//
// [[Rcpp::export]]
DataFrame insert_hinges(NumericVector w_length, NumericVector s_irrad, NumericVector hinges) {

    vector<double> w_length_in = Rcpp::as< vector<double> >(w_length);
    vector<double> s_irrad_in = Rcpp::as< vector<double> >(s_irrad);
    vector<double> hinges_in = Rcpp::as< vector<double> >(sort_unique(hinges));

    vector<double>::iterator h_low = lower_bound(hinges_in.begin(), hinges_in.end(), w_length_in.front());
    vector<double>::iterator h_high = lower_bound(hinges_in.begin(), hinges_in.end(), w_length_in.back());
#ifdef DEBUG
  cout << "h_low idx: " << distance(hinges_in.begin(), h_low) <<  "\n";
  cout << "h_high idx: " << distance(hinges_in.begin(), h_high) <<  "\n";
#endif

    vector<double>::iterator w_length_in_it = w_length_in.begin();
    vector<double>::iterator s_irrad_in_it = s_irrad_in.begin();

    vector<double>::size_type max_len = w_length_in.size() + distance(h_low, h_high);
    vector<double> w_length_out (max_len); // w_length_out.reserve(max_len);
    vector<double> s_irrad_out (max_len); //s_irrad_out.reserve(max_len);

    vector<double>::iterator w_length_out_it = w_length_out.begin();
    vector<double>::iterator s_irrad_out_it = s_irrad_out.begin();

    vector<double>::iterator wl_next_chunk_begin_it;
    vector<double>::iterator si_next_chunk_begin_it = s_irrad_in_it;

    for (vector<double>::iterator it=h_low; it<h_high; it++) {
      wl_next_chunk_begin_it = lower_bound(w_length_in_it-1, w_length_in.end(), *it);
      vector<double>::size_type num_elements = distance(w_length_in_it, wl_next_chunk_begin_it);
      si_next_chunk_begin_it = si_next_chunk_begin_it + num_elements;
#ifdef DEBUG
  cout << "num_elements: " << num_elements <<  "\n";
  cout << "it idx: " << distance(hinges_in.begin(), it) <<  "\n";
#endif
      w_length_out_it = copy(w_length_in_it, wl_next_chunk_begin_it, w_length_out_it);
      s_irrad_out_it = copy(s_irrad_in_it, si_next_chunk_begin_it, s_irrad_out_it);
#ifdef DEBUG
  cout << "w_length_out.size(): " << w_length_out.size() <<  "\n";
  cout << "s_irrad_out.size(): " << s_irrad_out.size() <<  "\n";
#endif
      if (*it < *wl_next_chunk_begin_it) {
        *w_length_out_it = *it;
         double interpol_irrad = *(si_next_chunk_begin_it-1) +
               (*si_next_chunk_begin_it - *(si_next_chunk_begin_it-1)) /
               (*wl_next_chunk_begin_it - *(wl_next_chunk_begin_it-1)) * (*it - *(wl_next_chunk_begin_it-1));
        *s_irrad_out_it = interpol_irrad;
        ++w_length_out_it;
        ++s_irrad_out_it;
      }
    w_length_in_it = wl_next_chunk_begin_it;
    s_irrad_in_it = si_next_chunk_begin_it;
    }
    w_length_out_it = copy(w_length_in_it, w_length_in.end(), w_length_out_it);
    s_irrad_out_it = copy(s_irrad_in_it, s_irrad_in.end(), s_irrad_out_it);
    w_length_out.erase(w_length_out_it,  w_length_out.end());
    s_irrad_out.erase(s_irrad_out_it,  s_irrad_out.end());
#ifdef DEBUG
  cout << "w_length_in.size(): " << w_length_in.size() <<  "\n";
  cout << "s_irrad_in.size(): " << s_irrad_in.size() <<  "\n";
  cout << "w_length_out.size(): " << w_length_out.size() <<  "\n";
  cout << "s_irrad_out.size(): " << s_irrad_out.size() <<  "\n";
#endif

    return DataFrame::create(_["w.length"] = w_length_out, _["s.irrad"] = s_irrad_out);
}
