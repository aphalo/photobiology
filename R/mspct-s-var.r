#' Variance of a collection of spectra
#'
#' Method to compute the "parallel" variance of values across members of a
#' collections of spectra or of a spectral object containing multiple spectra in
#' long form.
#'
#' @details Method specializations compute the variance at each wavelength
#'   across a group of spectra stored in an object of one of the classes defined
#'   in package 'photobiology'. Omission of NAs is done separately at each
#'   wavelength. Interpolation is not applied, so all spectra in \code{x} must
#'   share the same set of wavelengths. An error is triggered if this condition
#'   is nor fulfilled.
#'
#' @inheritParams s_mean
#'
#' @inherit s_mean note return
#'
#' @seealso See \code{\link[stats]{cor}} for details about \code{var()},
#'   which is used for the computations.
#'
#' @export
#'
#' @examples
#' s_var(sun_evening.mspct)
#'
s_var <- function(x, na.rm, ...)
  UseMethod("s_var")

#' @rdname s_var
#'
#' @export
#'
s_var.default <- function(x, na.rm = FALSE, ...) {
  warning("Metod 's_var()' not implementd for objects of class ",
          class(x)[1],
          ".")
  ifelse(is.any_mspct(x), generic_spct(), NA)
}

#' @rdname s_var
#'
#' @export
#'
s_var.generic_spct <- function(x, na.rm = FALSE, ...) {
  s_var(subset2mspct(x), na.rm = na.rm, ...)
}

#' @rdname s_var
#'
#' @export
#'
s_var.filter_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_filter(
    x = x,
    .fun = stats::var,
    na.rm = na.rm,
    col.name.tag = ".var",
    .fun.name = "Variance for"
  )
}

#' @rdname s_var
#'
#' @export
#'
s_var.source_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_source(
    x = x,
    .fun = stats::var,
    na.rm = na.rm,
    col.name.tag = ".var",
    .fun.name = "Variance for"
  )
}

#' @rdname s_var
#'
#' @export
#'
s_var.response_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_response(
    x = x,
    .fun = stats::var,
    na.rm = na.rm,
    col.name.tag = ".var",
    .fun.name = "Variance for"
  )
}

#' @rdname s_var
#'
#' @export
#'
s_var.reflector_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_reflector(
    x = x,
    .fun = stats::var,
    na.rm = na.rm,
    col.name.tag = ".var",
    .fun.name = "Variance for"
  )
}

#' @rdname s_var
#'
#' @export
#'
s_var.calibration_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_calibration(
    x = x,
    .fun = stats::var,
    na.rm = na.rm,
    col.name.tag = ".var",
    .fun.name = "Variance for"
  )
}

#' @rdname s_var
#'
#' @export
#'
s_var.cps_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_cps(
    x = x,
    .fun = stats::var,
    na.rm = na.rm,
    col.name.tag = ".var",
    .fun.name = "Variance for"
  )
}

#' @rdname s_var
#'
#' @export
#'
s_var.raw_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_raw(
    x = x,
    .fun = stats::var,
    na.rm = na.rm,
    col.name.tag = ".var",
    .fun.name = "Variance for"
  )
}
