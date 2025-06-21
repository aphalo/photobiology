#' Range of a collection of spectra
#'
#' Method to compute the "parallel" range of values across members of a
#' collection of spectra or of a spectral object containing multiple spectra in
#' long form.
#'
#' @details Method specializations compute the range at each wavelength across a
#'   group of spectra stored in an object of one of the classes defined in
#'   package 'photobiology'. Omission of NAs is done
#'   separately at each wavelength. Interpolation is not applied, so all spectra
#'   in \code{x} must share the same set of wavelengths. An error is triggered
#'   if this condition is nor fulfilled.
#'
#' @param x An R object. Currently this package defines methods for collections
#'    of spectral objects.
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
#' @seealso See \code{\link[base]{Extremes}} details on the \code{min()} and
#'   \code{max()} methods used for the computations.
#'
#' @export
#'
#' @examples
#' s_range(sun_evening.mspct)
#'
s_range <- function(x, na.rm, ...)
  UseMethod("s_range")

#' @describeIn s_range
#'
#' @export
#'
s_range.default <- function(x, na.rm = FALSE, ...) {
  warning("Metod 's_range()' not implementd for objects of class ",
          class(x)[1],
          ".")
  ifelse(is.any_mspct(x), generic_spct(), NA)
}

#' @describeIn s_range
#'
#' @export
#'
s_range.generic_spct <- function(x, na.rm = FALSE, ...) {
  s_range(subset2mspct(x), na.rm = na.rm, ...)
}

#' @describeIn s_range
#'
#' @export
#'
s_range.filter_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_filter(
    x,
    .fun = c(base::min, base::max),
    na.rm = na.rm,
    col.name.tag = c(".min", ".max"),
    .fun.name = "s_range of"
  )
}

#' @describeIn s_range
#'
#' @export
#'
s_range.source_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_source(
    x,
    .fun = c(base::min, base::max),
    na.rm = na.rm,
    col.name.tag = c(".min", ".max"),
    .fun.name = "s_range of"
  )
}

#' @describeIn s_range
#'
#' @export
#'
s_range.response_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_response(
    x,
    .fun = c(base::min, base::max),
    na.rm = na.rm,
    col.name.tag = c(".min", ".max"),
    .fun.name = "s_range of"
  )
}

#' @describeIn s_range
#'
#' @export
#'
s_range.reflector_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_reflector(
    x,
    .fun = c(base::min, base::max),
    na.rm = na.rm,
    col.name.tag = c(".min", ".max"),
    .fun.name = "s_range of"
  )
}

#' @describeIn s_range
#'
#' @export
#'
s_range.calibration_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_calibration(
    x,
    .fun = c(base::min, base::max),
    na.rm = na.rm,
    col.name.tag = c(".min", ".max"),
    .fun.name = "s_range of"
  )
}

#' @describeIn s_range
#'
#' @export
#'
s_range.cps_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_cps(
    x,
    .fun = c(base::min, base::max),
    na.rm = na.rm,
    col.name.tag = c(".min", ".max"),
    .fun.name = "s_range of"
  )
}

#' @describeIn s_range
#'
#' @export
#'
s_range.raw_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_raw(
    x,
    .fun = c(base::min, base::max),
    na.rm = na.rm,
    col.name.tag = c(".min", ".max"),
    .fun.name = "s_range of"
  )
}
