#' Median of a collection of spectra
#'
#' A method to compute the median of values across members of a collections of
#' spectra. Computes the median at each wavelength across all the spectra in the
#' collection returning a spectral object.
#'
#' @param x An R object. Currently this package defines methods for collections
#'   of spectral objects.
#' @param na.rm	logical. A value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param ...	Further arguments passed to or from other methods.
#'
#' @return If \code{x} is a collection spectral of objects, such as a
#'   "filter_mspct" object, the returned object is of same class as the members
#'   of the collection, such as "filter_spct", containing the median spectrum.
#'
#' @note Omission of NAs is done separately at each wavelength. Interpolation is
#'   not applied, so all spectra in \code{x} must share the same set of
#'   wavelengths.
#'
#'   Objects of classes raw_spct and cps_spct can contain data from multiple
#'   scans. This functions are implemented for these classes only for the case
#'   when all member spectra contain data for a single scan, or spliced into a
#'   single column in the case of cps_spct members.
#'
#' @seealso See \code{\link[stats]{median}} for the \code{median()} method used
#'   for the computations.
#'
#' @export
#'
s_median <- function(x, na.rm, ...) UseMethod("s_median")

#' @describeIn s_median
#'
#' @export
#'
s_median.default <- function(x, na.rm = FALSE, ...) {
  warning("Metod 's_median()' not implementd for objects of class ", class(x)[1], ".")
  ifelse(is.any_mspct(x), do.call(class(x[[1]])[1], args = list()), NA)
}

#' @describeIn s_median
#'
#' @export
#'
s_median.generic_spct <- function(x, na.rm = FALSE, ...) {
  if (getMultipleWl(x) > 1) {
    s_median(subset2mspct(x), na.rm = na.rm, ...)
  } else {
    x
  }
}

#' @describeIn s_median
#'
#' @export
#'
s_median.source_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_source(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}

#' @describeIn s_median
#'
#' @export
#'
s_median.response_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_response(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}

#' @describeIn s_median
#'
#' @export
#'
s_median.filter_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_filter(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}

#' @describeIn s_median
#'
#' @export
#'
s_median.reflector_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_reflector(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}

#' @describeIn s_median
#'
#' @export
#'
s_median.calibration_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_calibration(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}

#' @describeIn s_median
#'
#' @export
#'
s_median.cps_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_cps(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}

#' @describeIn s_median
#'
#' @export
#'
s_median.raw_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_raw(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}
