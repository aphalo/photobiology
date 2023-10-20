#' Divide two spectra, even if the wavelengths values differ
#'
#' The wavelength vectors of the two spectra are merged, and the missing
#' spectral values are calculated by interpolation. After this, the two spectral
#' values at each wavelength are operated upon.
#'
#' @param w.length1 numeric vector of wavelength (nm) of denominator.
#' @param w.length2 numeric vector of wavelength (nm) of divisor.
#' @param s.irrad1 a numeric vector of spectral values of denominator.
#' @param s.irrad2 a numeric vector of spectral values of divisor.
#' @param trim a character string with value "union" or "intersection".
#' @param na.rm a logical value, if TRUE, not the default, NAs in the input are
#'   replaced with zeros.
#'
#' @return a dataframe with two numeric variables.
#' \item{w.length}{A numeric vector with the wavelengths (nm) obtained
#'   by "fusing" w.length1 and
#'   w.length2. w.length contains all the unique vales, sorted in ascending
#'   order.}
#'   \item{s.irrad}{A numeric vector with the ratio between the two spectral
#'   values at each wavelength.}
#'
#' @details If trim=="union" spectral values are calculated for the whole range
#'   of wavelengths covered by at least one of the input spectra, and missing
#'   values are set in each input spectrum to zero before addition. If
#'   trim=="intersection" then the range of wavelengths covered by both input
#'   spectra is returned, and the non-overlapping regions discarded. If
#'   w.length2==NULL, it is assumed that both spectra are measured at the same
#'   wavelengths, and a simple addition is used, ensuring fast calculation.
#' @export
#'
#' @examples
#'
#' head(sun.data)
#' one.data <- with(sun.data, div_spectra(w.length, w.length, s.e.irrad, s.e.irrad))
#' head(one.data)
#' tail(one.data)
#'
#' @family low-level functions operating on numeric vectors.
#'
div_spectra <-
  function(w.length1,
           w.length2 = NULL,
           s.irrad1,
           s.irrad2,
           trim = "union",
           na.rm = FALSE) {
    return(
      oper_spectra(
        w.length1 = w.length1,
        w.length2 = w.length2,
        s.irrad1 = s.irrad1,
        s.irrad2 = s.irrad2,
        trim = "union",
        na.rm = FALSE,
        bin.oper = `/`
      )
    )
  }
