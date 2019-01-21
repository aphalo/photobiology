#' Standard Deviation of a collection of spectra
#'
#' A method to compute the standard deviation of values across members of a
#' collections of spectra. Computes the standard deviation at each wavelength
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
#'   containing the standard deviation among the spectra at each wavelength
#'   in a column with name ending in ".sd".
#'
#' @note Omission of NAs is done separately at each wavelength. Interpolation is
#'   not applied, so all spectra in \code{x} must share the same set of
#'   wavelengths.
#'
#' @seealso See \code{\link[stats]{sd}} for details about \code{sd()} methods
#'   for other classes.
#'
#' @export
#'
s_sd <- function(x, na.rm, ...) UseMethod("s_sd")

#' @describeIn s_sd
#'
#' @export
#'
s_sd.default <- function(x, na.rm = FALSE, ...) {
  warning("Metod 'sd()' not implementd for objects of class ", class(x)[1], ".")
  ifelse(is.any_mspct(x), generic_spct(), NA)
}

#' @describeIn s_sd
#'
#' @export
#'
s_sd.filter_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_filter(x = x, .fun = stats::sd, na.rm = na.rm, col.name.tag = ".sd", .fun.name = "Standard deviation for")
}

#' @describeIn s_sd
#'
#' @export
#'
s_sd.source_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_source(x = x, .fun = stats::sd, na.rm = na.rm, col.name.tag = ".sd", .fun.name = "Standard deviation for")
}

#' @describeIn s_sd
#'
#' @export
#'
s_sd.response_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_response(x = x, .fun = stats::sd, na.rm = na.rm, col.name.tag = ".sd", .fun.name = "Standard deviation for")
}

#' @describeIn s_sd
#'
#' @export
#'
s_sd.reflector_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_reflector(x = x, .fun = stats::sd, na.rm = na.rm, col.name.tag = ".sd", .fun.name = "Standard deviation for")
}
