#' Median of a collection of spectra
#'
#' Method to compute the "parallel" median of values across members of a
#' collection of spectra or of a spectral object containing multiple spectra in
#' long form.
#'
#' @details Method specializations compute the median at each wavelength across
#'   a group of spectra stored in an object of one of the classes defined in
#'   package 'photobiology'. Omission of NAs is done
#'   separately at each wavelength. Interpolation is not applied, so all spectra
#'   in \code{x} must share the same set of wavelengths. An error is triggered
#'   if this condition is nor fulfilled.
#'
#' @inheritParams s_mean
#'
#' @inherit s_mean note return
#'
#' @inheritSection s_mean Deepest Curves
#'
#' @seealso See \code{\link[stats]{median}} for the \code{median()} method used
#'   for the computations.
#'
#' @export
#'
#' @examples
#' s_median(sun_evening.mspct)
#'
s_median <- function(x, na.rm, ...) UseMethod("s_median")

#' @rdname s_median
#'
#' @export
#'
s_median.default <- function(x, na.rm = FALSE, ...) {
  warning("Metod 's_median()' not implementd for objects of class ", class(x)[1], ".")
  ifelse(is.any_mspct(x), do.call(class(x[[1]])[1], args = list()), NA)
}

#' @rdname s_median
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

#' @rdname s_median
#'
#' @export
#'
s_median.source_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_source(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}

#' @rdname s_median
#'
#' @export
#'
s_median.response_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_response(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}

#' @rdname s_median
#'
#' @export
#'
s_median.filter_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_filter(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}

#' @rdname s_median
#'
#' @export
#'
s_median.reflector_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_reflector(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}

#' @rdname s_median
#'
#' @export
#'
s_median.calibration_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_calibration(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}

#' @rdname s_median
#'
#' @export
#'
s_median.cps_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_cps(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}

#' @rdname s_median
#'
#' @export
#'
s_median.raw_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_raw(x = x, .fun = stats::median, na.rm = na.rm, .fun.name = "Median of")
}
