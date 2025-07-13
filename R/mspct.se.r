#' Standard Error of a collection of spectra
#'
#' Method to compute the "parallel" standard error of the mean across members of
#' a collection of spectra or of a spectral object containing multiple spectra
#' in long form.
#'
#' @details Method specializations compute the standard error of the mean at
#'   each wavelength across a group of spectra stored in an object of one of the
#'   classes defined in package 'photobiology'. Omission of NAs is done
#'   separately at each wavelength. Interpolation is not applied, so all spectra
#'   in \code{x} must share the same set of wavelengths. An error is triggered
#'   if this condition is nor fulfilled.
#'
#' @inheritParams s_mean
#'
#' @inherit s_mean note return
#'
#' @seealso See \code{\link[stats]{sd}} for details about \code{sd()} methods
#'   for other classes.
#'
#' @export
#'
#' @examples
#' s_se(sun_evening.mspct)
#'
s_se <- function(x, na.rm, ...)
  UseMethod("s_se")

#' @rdname s_se
#'
#' @export
#'
s_se.default <- function(x, na.rm = FALSE, ...) {
  warning("Metod 'se()' not implementd for objects of class ",
          class(x)[1],
          ".")
  ifelse(is.any_mspct(x), generic_spct(), NA)
}

#' @rdname s_se
#'
#' @export
#'
s_se.generic_spct <- function(x, na.rm = FALSE, ...) {
  s_se(subset2mspct(x), na.rm = na.rm, ...)
}

#' @rdname s_se
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

#' @rdname s_se
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

#' @rdname s_se
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

#' @rdname s_se
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

#' @rdname s_se
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

#' @rdname s_se
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

#' @rdname s_se
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
