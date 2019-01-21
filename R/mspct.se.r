#' Standard Error of a collection of spectra
#'
#' A method to compute the standard error of values across members of a
#' collections of spectra. Computes the standard error at each wavelength
#' across all the spectra in the collection returning a spectral object.
#'
#' @param x An R object. Currently this package defines methods for collections
#'   of spectral objects.
#' @param na.rm	logical. A value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param ...	Further arguments passed to or from other methods.
#'
#' @return If \code{x} is a collection spectral of objects, such as a
#'   "filter_mspct" object, the returned object is of class "generic_spct",
#'   containing the standard error among the spectra at each wavelength
#'   in a column with name ending in ".se".
#'
#' @note Omission of NAs is done separately at each wavelength. Interpolation is
#'   not applied, so all spectra in \code{x} must share the same set of
#'   wavelengths.
#'
#' @export
#'
s_se <- function(x, na.rm, ...) UseMethod("s_se")

#' @describeIn s_se
#'
#' @export
#'
s_se.default <- function(x, na.rm = FALSE, ...) {
  warning("Metod 'se()' not implementd for objects of class ", class(x)[1], ".")
  ifelse(is.any_mspct(x), generic_spct(), NA)
}

#' @describeIn s_se
#'
#' @export
#'
s_se.filter_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_filter(x = x, .fun = se, na.rm = na.rm, col.name.tag = ".se", .fun.name = "Standard error for")
}

#' @describeIn s_se
#'
#' @export
#'
s_se.source_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_source(x = x, .fun = se, na.rm = na.rm, col.name.tag = ".se", .fun.name = "Standard error for")
}

#' @describeIn s_se
#'
#' @export
#'
s_se.response_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_response(x = x, .fun = se, na.rm = na.rm, col.name.tag = ".se", .fun.name = "Standard error for")
}

#' @describeIn s_se
#'
#' @export
#'
s_se.reflector_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_reflector(x = x, .fun = se, na.rm = na.rm, col.name.tag = ".se", .fun.name = "Standard error for")
}
