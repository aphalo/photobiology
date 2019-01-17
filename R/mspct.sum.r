#' Sum from collection of spectra
#'
#' A method to compute the sum of values across members of a collections of
#' spectra. Computes the sum at each wavelength across all the spectra in the
#' collection returning a spectral object.
#'
#' @param x An R object. Currently this package defines methods for collections of
#'    spectral objects.
#' @param na.rm	logical. A value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param ...	Further arguments passed to or from other methods.
#'
#' @return If \code{x} is a collection spectral of objects, such as a
#'   "filter_mspct" object, the returned object is of same class as the
#'   members of the collection, such as "filter_spct", containing the sum
#'   of the spectra.
#'
#' @note Omission of NAs is done separately at each wavelength. Interpolation is
#'   not applied, so all spectra in \code{x} must share the same set of
#'   wavelengths.
#'
#' A sum of transmitances or reflectances is no longer a well defined
#' physical quanttiy, and these sum operations return an object of class
#' generic_spct.
#'
#' @seealso See \code{\link[base]{sum}} for the \code{sum()} method used for
#'   the computations.
#'
#' @export
#'
s_sum <- function(x, na.rm, ...) UseMethod("s_sum")

#' @describeIn s_sum
#'
#' @export
#'
s_sum.default <- function(x, na.rm = FALSE, ...) {
  warning("Metod 's_sum()' not implementd for objects of class ", class(x)[1], ".")
  ifelse(is.any_mspct(x), do.call(class(x[[1]])[1], args = list()), NA)
}

#' @describeIn s_sum
#'
#' @export
#'
s_sum.filter_mspct <- function(x, na.rm = FALSE, ...) {
  warning("A sum of Tfr values does not yield Tfr values, while a summ of A values yields A values!!")
  rowwise_filter(x, .fun = base::sum, na.rm = na.rm, col.name.tag = ".sum", .fun.name = "Sum of")
}

#' @describeIn s_sum
#'
#' @export
#'
s_sum.source_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_source(x, .fun = base::sum, na.rm = na.rm, .fun.name = "Sum of")
}

#' @describeIn s_sum
#'
#' @export
#'
s_sum.response_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_response(x, .fun = base::sum, na.rm = na.rm, .fun.name = "Sum of")
}

#' @describeIn s_sum
#'
#' @export
#'
s_sum.reflector_mspct <- function(x, na.rm = FALSE, ...) {
  warning("A sum of Rfr values does not yield Rfr values!!")
  rowwise_reflector(x, .fun = base::sum, na.rm = na.rm, col.name.tag = ".sum", .fun.name = "Sum of")
}

#' @describeIn s_sum
#'
#' @export
#'
s_sum.calibration_mspct <- function(x, na.rm = FALSE, ...) {
  warning("A sum of irrad.mult values does not yield irrad.mult values!!")
  rowwise_calibration(x, .fun = base::sum, na.rm = na.rm, col.name.tag = ".sum", .fun.name = "Sum of")
}
