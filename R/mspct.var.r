#' Variance of a collection of spectra
#'
#' A method to compute the variance of values across members of a collections of
#' spectra. Computes the variance at each wavelength across all the spectra in
#' the collection returning a spectral object.
#'
#' Variance method for collections of spectra. Computes the variance at each
#' wavelength across all the spectra in the collection.
#'
#' @param x An R object. Currently this package defines methods for collections
#'   of spectral objects.
#' @param na.rm	logical. A value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param ...	Further arguments passed to or from other methods.
#'
#' @return If \code{x} is a collection spectral of objects, such as a
#'   "filter_mspct" object, the returned object is of class "generic_spct",
#'   containing the variance among the spectra at each wavelength
#'   in a column with name ending in ".var".
#'
#' @note Omission of NAs is done separately at each wavelength. Interpolation is
#'   not applied, so all spectra in \code{x} must share the same set of
#'   wavelengths.
#'
#' @seealso See \code{\link[stats]{cor}} for details about \code{var()}, which
#'   is used for the computations.
#'
#' @export
#'
s_var <- function(x, na.rm, ...) UseMethod("s_var")

#' @describeIn s_var
#'
#' @export
#'
s_var.default <- function(x, na.rm = FALSE, ...) {
  warning("Metod 's_var()' not implementd for objects of class ", class(x)[1], ".")
  ifelse(is.any_mspct(x), generic_spct(), NA)
}

#' @describeIn s_var
#'
#' @export
#'
s_var.filter_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_filter(x = x, .fun = stats::var, na.rm = na.rm, col.name.tag = ".var", .fun.name = "Variance for")
}

#' @describeIn s_var
#'
#' @export
#'
s_var.source_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_source(x = x, .fun = stats::var, na.rm = na.rm, col.name.tag = ".var", .fun.name = "Variance for")
}

#' @describeIn s_var
#'
#' @export
#'
s_var.response_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_response(x = x, .fun = stats::var, na.rm = na.rm, col.name.tag = ".var", .fun.name = "Variance for")
}

#' @describeIn s_var
#'
#' @export
#'
s_var.reflector_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_reflector(x = x, .fun = stats::var, na.rm = na.rm, col.name.tag = ".var", .fun.name = "Variance for")
}

#' @describeIn s_var
#'
#' @export
#'
s_var.calibration_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_calibration(x = x, .fun = stats::var, na.rm = na.rm, col.name.tag = ".var", .fun.name = "Variance for")
}
