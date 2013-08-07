#' Sanity check of a spectrum.
#'
#' This function checks a radiation spectrum for compliance with assumptions used in calculations.
#'
#' @usage check_spectrum(w.length, s.irrad)
#' 
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral (energy) irradiances (W m-2 nm-1)
#' 
#' @return a single logicalk value indicating whether test was passed or not
#' @export
#' 
#' @keywords manip misc
#'
#' @examples
#' data(sun.data)
#' with(sun.data, check_spectrum(w.length, s.e.irrad))
#' 
check_spectrum <- function(w.length, s.irrad) {
  # check for NAs
  if (any(is.na(w.length)|is.na(s.irrad))){
    warning("Error: at least one NA value in wavelengths vector and/or s.e.irrad vector")
    return(FALSE)
  }
  if (is.unsorted(w.length, strictly=TRUE)) {
    warning("Error: wavelengths should be sorted in ascending order")
    return(FALSE)
  }
  if (length(w.length) != length(s.irrad)){
    warning("Error: wavelengths vector and s.e.irrad vector should have same length")
    return(FALSE)
  }
  # warn if w.length values are not reasonable
  if (min(w.length < 200.0) || max(w.length > 1000.0)){
    warning("Warning: wavelength values should be in nm\n data contains values < 200 nm and/or > 1000 nm")
  }
  # test average wavelength delta
  w.length.resolution <- (max(w.length) - min(w.length))/ length(w.length)
  if (w.length.resolution > 10.0) {
    warning("Average wavelength resolution > 10 nm, too low for accurate results.")
    return(FALSE)
  } else if (w.length.resolution > 5.0) {
    warning("Average wavelength resolution > 5 nm, calculations may be inaccurate in some cases.")
  } 
  return(TRUE)
}

