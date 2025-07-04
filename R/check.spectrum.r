#' Sanity check a spectrum
#'
#' Checks spectral irradiance data in \code{numeric} vectors for
#' compliance with assumptions used in calculations.
#'
#' @param w.length numeric vector of wavelengths [\eqn{nm}].
#' @param s.irrad numeric Corresponding vector of spectral (energy) irradiances
#'   [\eqn{W\,m^{-2}\,nm^{-1}}{W m-2 nm-1}].
#'
#' @return A single \code{logical} value indicating whether test was passed or
#'   not
#' @export
#'
#' @examples
#' with(sun.data, check_spectrum(w.length, s.e.irrad))
#'
#' @family data validity check functions
#'
check_spectrum <- function(w.length, s.irrad) {
  pass <- TRUE
  # check for NAs
  if (any(is.na(w.length) | is.na(s.irrad))) {
    warning("Error: at least one NA value in wavelengths vector ",
            "and/or s.e.irrad vector")
    pass <- FALSE
  }
  if (is.unsorted(w.length, na.rm = TRUE, strictly = TRUE)) {
    warning("Error: wavelengths should be sorted in ascending order")
    pass <- FALSE
  }
  if (length(w.length) != length(s.irrad)) {
    warning("Error: wavelengths vector and s.e.irrad vector ",
            "should have same length")
    pass <- FALSE
  }
  # warn if w.length values are not reasonable
  if (min(w.length < 100.0) || max(w.length > 6000.0)) {
    warning("Wavelength values should be in nm. ",
            "Spectrum contains values < 100 nm and/or > 6000 nm")
  }
  # test average wavelength delta
  w.length.resolution <- (max(w.length) - min(w.length))/ length(w.length)
  if (w.length.resolution > 10.0) {
    warning("Average wavelength resolution > 10 nm, ",
            "too low for accurate results.")
    pass <- FALSE
  } else if (w.length.resolution > 5.0) {
    warning("Average wavelength resolution > 5 nm, calculations ",
            "can be inaccurate with steep slopes, narrow peaks or valleys.")
  }
  return(pass)
}

#' Sanity check of wavelengths (internal function).
#'
#' This function checks a w.length vector for compliance with assumptions
#' expected for valid calculations.
#'
#' @param w.length numeric array of wavelength (nm)
#'
#' @return a single logical value indicating whether test was passed or not
#'
#' @examples
#'
#' with(sun.data, photobiology:::check_w.length(w.length))
#'
#' @family data validity check functions
#'
check_w.length <- function(w.length) {
  pass <- TRUE
  if (is.unsorted(w.length, strictly=TRUE)) {
    warning("Possible error: wavelengths are not sorted in ascending order")
    pass <- FALSE
  }
  # warn if w.length values are not reasonable
  if (min(w.length < 100.0) || max(w.length > 6000.0)){
    warning("Wavelength values should be in nm. ",
            "Spectrum contains values < 100 nm and/or > 6000 nm")
    pass <- FALSE
  }
  # test average wavelength delta
  w.length.resolution <- (max(w.length) - min(w.length))/ length(w.length)
  if (w.length.resolution > 10.0) {
     warning("Average wavelength resolution > 10 nm, ",
            "too low for accurate results.")
    return(FALSE)
    pass <- FALSE
  } else if (w.length.resolution > 5.0) {
     warning("Average wavelength resolution > 5 nm, calculations ",
            "can be inaccurate with steep slopes, narrow peaks or valleys.")
  }
  return(pass)
}
