#' Range of a collection of spectra
#'
#' A method to compute the range of values across members of a collections of
#' spectra. Computes the max and min at each wavelength across all the spectra
#' in the collection returning a spectral object.
#'
#' @param x An R object. Currently this package defines methods for collections of
#'    spectral objects.
#' @param na.rm	logical. A value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param ...	Further arguments passed to or from other methods.
#'
#' @return If \code{x} is a collection spectral of objects, such as a
#'   "filter_mspct" object, the returned object is of same class as the
#'   members of the collection, such as "filter_spct", containing the mean
#'   spectrum.
#'
#' @note Trimming of extreme values and omission of NAs is done separately at
#'   each wavelength. Interpolation is not applied, so all spectra in \code{x}
#'   must share the same set of wavelengths.
#'
#'   Objects of classes raw_spct and cps_spct can contain data from multiple
#'   scans. This functions are implemented for these classes only for the case
#'   when all member spectra contain data for a single scan, or spliced into a
#'   single column in the case of cps_spct members.
#'
#' @seealso See \code{\link[base]{Extremes}} details on the \code{min()} and
#'   \code{max()} methods used for the computations.
#'
#' @export
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
