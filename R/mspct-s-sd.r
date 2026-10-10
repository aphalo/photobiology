#' Standard Deviation of a collection of spectra
#'
#' Method to compute the "parallel" standard deviation of values across members
#' of a collection of spectra or of a spectral object containing multiple
#' spectra in long form.
#'
#' @details Method specializations compute the standard deviation at each
#'   wavelength across a group of spectra stored in an object of one of the
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
#' s_sd(sun_evening.mspct)
#'
s_sd <- function(x, na.rm, ...)
  UseMethod("s_sd")

#' @rdname s_sd
#'
#' @export
#'
s_sd.default <- function(x, na.rm = FALSE, ...) {
  warning("Metod 'sd()' not implementd for objects of class ",
          class(x)[1],
          ".")
  ifelse(is.any_mspct(x), generic_spct(), NA)
}

#' @rdname s_sd
#'
#' @export
#'
s_sd.generic_spct <- function(x, na.rm = FALSE, ...) {
  s_sd(subset2mspct(x), na.rm = na.rm, ...)
}

#' @rdname s_sd
#'
#' @export
#'
s_sd.filter_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_filter(
    x = x,
    .fun = stats::sd,
    na.rm = na.rm,
    col.name.tag = ".sd",
    .fun.name = "Standard deviation for"
  )
}

#' @rdname s_sd
#'
#' @export
#'
s_sd.source_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_source(
    x = x,
    .fun = stats::sd,
    na.rm = na.rm,
    col.name.tag = ".sd",
    .fun.name = "Standard deviation for"
  )
}

#' @rdname s_sd
#'
#' @export
#'
s_sd.response_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_response(
    x = x,
    .fun = stats::sd,
    na.rm = na.rm,
    col.name.tag = ".sd",
    .fun.name = "Standard deviation for"
  )
}

#' @rdname s_sd
#'
#' @export
#'
s_sd.reflector_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_reflector(
    x = x,
    .fun = stats::sd,
    na.rm = na.rm,
    col.name.tag = ".sd",
    .fun.name = "Standard deviation for"
  )
}

#' @rdname s_sd
#'
#' @export
#'
s_sd.calibration_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_calibration(
    x = x,
    .fun = stats::sd,
    na.rm = na.rm,
    col.name.tag = ".sd",
    .fun.name = "Standard deviation for"
  )
}

#' @rdname s_sd
#'
#' @export
#'
s_sd.cps_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_cps(
    x = x,
    .fun = stats::sd,
    na.rm = na.rm,
    col.name.tag = ".sd",
    .fun.name = "Standard deviation for"
  )
}

#' @rdname s_sd
#'
#' @export
#'
s_sd.raw_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_raw(
    x = x,
    .fun = stats::sd,
    na.rm = na.rm,
    col.name.tag = ".sd",
    .fun.name = "Standard deviation for"
  )
}
