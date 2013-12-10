#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// #define DEBUG
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
//' @usage insert_hingesC(w_length, s_irrad, hinges) 
//' 
//' @param w_length numeric array of wavelength (nm)
//' @param s_irrad numeric array of spectral irradiance values
//' @param hinges a numeric array giving the wavelengths at which the s.irrad should be inserted by
//' interpolation, no interpolation is indicated by an empty array (numeric(0))
//' 
//' @return a data.frame with variables \code{w.length} and \code{s.irrad}
//' @name insert_hingesC
//' @keywords manip misc
//' @export
//' @examples
//' with(sun.data, insert_hingesC(w.length, s.e.irrad, c(399.99,400.00,699.99,700.00)))
//' with(sun.data, insert_hingesC(w.length, s.e.irrad, c(199.99,200.00,399.50,399.99,400.00,699.99,700.00,799.99,1000.00)))
//
// I have tried to optimize this function as much as possible.
// It assumes that hinges are sorted in increasing order.
// It uses binary search to find the insertion point, and also does this search in the w_length
// vector only in the "tail" starting at the previously inserted hinge.
// The second example tests several boundary conditions: response for hinges outside the range of the spectrum, two succesive
// insertions between the a single pair of values in the spectral data, and for insertion of
// hinges at walengths already present in the data. A test for the boundary condition
// when an empty numeric vector is passed to hinges fails. 
// with(sun.data, insert_hingesC(w.length, s.e.irrad, c())) in Rcpp
//
// [[Rcpp::export]]
DataFrame insert_hingesC (NumericVector w_length, NumericVector s_irrad, NumericVector hinges) {

    NumericVector::iterator h_low = lower_bound(hinges.begin(), hinges.end(), *w_length.begin());
    NumericVector::iterator h_high = lower_bound(hinges.begin(), hinges.end(), *(w_length.end() - 1));
    NumericVector::iterator inserted;
    NumericVector::iterator after_inserted;

      #ifdef DEBUG    
        cout << "h_low: " << distance(hinges.begin(), h_low) <<  "\n";
        cout << "h_high: " << distance(hinges.begin(), h_high) <<  "\n";
      #endif
    int i = 0;
    for (NumericVector::iterator it=h_low; it<h_high; it++) {
      after_inserted = lower_bound(w_length.begin() + i, w_length.end(), *it);
      #ifdef DEBUG    
        cout << "*after_inserted: " << *after_inserted << "; *it: " << *it << "\n";
        cout << "after_inserted: " << distance(w_length.begin(), after_inserted) << "\n";
      #endif

      if (*it < *after_inserted) {
        inserted = w_length.insert(after_inserted, *it);
      #ifdef DEBUG    
        cout << "inserted: " << distance(w_length.begin(), inserted) <<  "\n";
        cout << "*inserted: " << *inserted <<  "\n";
        cout << "*it: " << *it << "\n";
      #endif
        i = distance(w_length.begin(), inserted);
        double interpol_irrad = s_irrad[i-1] + (s_irrad[i] - s_irrad[i-1]) /  (w_length[i+1] - w_length[i-1]) * (w_length[i] - w_length[i-1]);
      #ifdef DEBUG    
        cout << "i: " << i <<  "\n";
        cout << "interpol_irrad: " << interpol_irrad <<  "\n";
        cout << "*it: " << *it << "\n";
      #endif
        s_irrad.insert(s_irrad.begin() + i, interpol_irrad);
      }
    }
    return DataFrame::create(Named("w.length") = w_length,
                             Named("s.irrad") = s_irrad);
}
