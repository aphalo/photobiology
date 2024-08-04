#' Standard Error of a collection of spectra
#'
#' Method to compute the "parallel" standard error of the mean across members of a
#' collection of spectra or of a spectral object containing multiple spectra in
#' long form.
#'
#' @details Method specializations compute the standard error of the mean at
#'   each wavelength across a group of spectra stored in an object of one of the
#'   classes defined in package 'photobiology'. Omission of NAs is done
#'   separately at each wavelength. Interpolation is not applied, so all spectra
#'   in \code{x} must share the same set of wavelengths. An error is triggered
#'   if this condition is nor fulfilled.
#'
#' @param x An R object. Currently this package defines methods for collections
#'   of spectral objects.
#' @param na.rm	logical. A value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param ...	Further arguments passed to or from other methods.
#'
#' @return If \code{x} is a collection spectral of objects, such as a
#'   \code{"filter_mspct"} object, the returned object is of same class as the
#'   members of the collection, such as \code{"filter_spct"}, containing the
#'   summary spectrum, with variables with names tagged for summaries other
#'   than mean or median.
#'
#' @note Objects of classes \code{raw_spct} and \code{cps_spct} can contain data
#'   from multiple scans in multiple variables or "columns". The methods accept
#'   as arguments objects of these classes only if spectra contain data for a
#'   single spectrometer scan. In the case of \code{cps_spct} objects, a single
#'   column can also contain data from multiple scans spliced into a single
#'   variable.
#'
#' @export
#'
#' @examples
#' s_se(sun_evening.mspct)
#'
s_se <- function(x, na.rm, ...)
  UseMethod("s_se")

#' @describeIn s_se
#'
#' @export
#'
s_se.default <- function(x, na.rm = FALSE, ...) {
  warning("Metod 'se()' not implementd for objects of class ",
          class(x)[1],
          ".")
  ifelse(is.any_mspct(x), generic_spct(), NA)
}

#' @describeIn s_se
#'
#' @export
#'
s_se.generic_spct <- function(x, na.rm = FALSE, ...) {
  s_se(subset2mspct(x), na.rm = na.rm, ...)
}

#' @describeIn s_se
#'
#' @export
#'
s_se.source_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_source(
    x = x,
    .fun = se,
    na.rm = na.rm,
    col.name.tag = ".se",
    .fun.name = "Standard error for"
  )
}

#' @describeIn s_se
#'
#' @export
#'
s_se.response_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_response(
    x = x,
    .fun = se,
    na.rm = na.rm,
    col.name.tag = ".se",
    .fun.name = "Standard error for"
  )
}

#' @describeIn s_se
#'
#' @export
#'
s_se.filter_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_filter(
    x = x,
    .fun = se,
    na.rm = na.rm,
    col.name.tag = ".se",
    .fun.name = "Standard error for"
  )
}

#' @describeIn s_se
#'
#' @export
#'
s_se.reflector_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_reflector(
    x = x,
    .fun = se,
    na.rm = na.rm,
    col.name.tag = ".se",
    .fun.name = "Standard error for"
  )
}

#' @describeIn s_se
#'
#' @export
#'
s_se.calibration_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_reflector(
    x = x,
    .fun = se,
    na.rm = na.rm,
    col.name.tag = ".se",
    .fun.name = "Standard error for"
  )
}

#' @describeIn s_se
#'
#' @export
#'
s_se.cps_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_cps(
    x = x,
    .fun = se,
    na.rm = na.rm,
    col.name.tag = ".se",
    .fun.name = "Standard error for"
  )
}

#' @describeIn s_se
#'
#' @export
#'
s_se.raw_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_raw(
    x = x,
    .fun = se,
    na.rm = na.rm,
    col.name.tag = ".se",
    .fun.name = "Standard error for"
  )
}
