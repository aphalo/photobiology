#' Mean and standard error from collection of spectra
#'
#' A method to compute the mean of values across members of a collections of
#' spectra. Computes the mean at each wavelength across all the spectra in the
#' collection returning a spectral object.
#'
#' @param x An R object Currently this package defines methods for collections of
#'    spectral objects.
#' @param na.rm	logical A value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param mult	numeric number of multiples of standard error.
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
#' @seealso See \code{\link[base]{mean}} for the \code{mean()} method used for
#'   the computations.
#'
#' @export
#'
s_mean_se <- function(x, na.rm, mult, ...)
  UseMethod("s_mean_se")

#' @describeIn s_mean_se
#'
#' @export
#'
s_mean_se.default <- function(x, na.rm = FALSE, mult = 1, ...) {
  warning("Metod 's_mean_se()' not implementd for objects of class ",
          class(x)[1],
          ".")
  ifelse(is.any_mspct(x), do.call(class(x[[1]])[1], args = list()), NA)
}

#' @describeIn s_mean_se
#'
#' @export
#'
s_mean_se.filter_mspct <-
  function(x, na.rm = FALSE, mult = 1, ...) {
    rowwise_filter(
      x,
      .fun = list(base::mean, se),
      col.name.tag = c("", ".se"),
      na.rm = na.rm,
      mult = mult,
      .fun.name = "Mean and SEM of"
    )
  }

#' @describeIn s_mean_se
#'
#' @export
#'
s_mean_se.source_mspct <-
  function(x, na.rm = FALSE, mult = 1, ...) {
    rowwise_source(
      x,
      .fun = list(base::mean, se),
      col.name.tag = c("", ".se"),
      na.rm = na.rm,
      mult = mult,
      .fun.name = "Mean and SEM of"
    )
  }

#' @describeIn s_mean_se
#'
#' @export
#'
s_mean_se.response_mspct <-
  function(x, na.rm = FALSE, mult = 1, ...) {
    rowwise_response(
      x,
      .fun = list(base::mean, se),
      col.name.tag = c("", ".se"),
      na.rm = na.rm,
      mult = mult,
      .fun.name = "Mean and SEM of"
    )
  }

#' @describeIn s_mean_se
#'
#' @export
#'
s_mean_se.reflector_mspct <-
  function(x, na.rm = FALSE, mult = 1, ...) {
    rowwise_reflector(
      x,
      .fun = list(base::mean, se),
      col.name.tag = c("", ".se"),
      na.rm = na.rm,
      mult = mult,
      .fun.name = "Mean and SEM of"
    )
  }

#' @describeIn s_mean_se
#'
#' @export
#'
s_mean_se.calibration_mspct <-
  function(x, na.rm = FALSE, mult = 1, ...) {
    rowwise_calibration(
      x,
      .fun = list(base::mean, se),
      col.name.tag = c("", ".se"),
      na.rm = na.rm,
      mult = mult,
      .fun.name = "Mean and SEM of"
    )
  }

#' @describeIn s_mean_se
#'
#' @export
#'
s_mean_se.cps_mspct <- function(x, na.rm = FALSE, mult = 1, ...) {
  rowwise_cps(
    x,
    .fun = list(base::mean, se),
    col.name.tag = c("", ".se"),
    na.rm = na.rm,
    mult = mult,
    .fun.name = "Mean and SEM of"
  )
}

#' @describeIn s_mean_se
#'
#' @export
#'
s_mean_se.raw_mspct <- function(x, na.rm = FALSE, mult = 1, ...) {
  rowwise_raw(
    x,
    .fun = list(base::mean, se),
    col.name.tag = c("", ".se"),
    na.rm = na.rm,
    mult = mult,
    .fun.name = "Mean and SEM of"
  )
}

# Helper function, not exported
#
#' Standard error of the mean
#'
#' @param x	numeric vector
#' @param mult	numeric number of multiples of standard error
#'
#' @note mult can be used to construct confidence intervals
#'
#' @keywords internal
#'
se <- function(x, na.rm = FALSE, mult = 1, ...) {
  stopifnot(is.numeric(x))
  if (na.rm) {
    x <- stats::na.omit(x)
  }
  mult * sqrt(stats::var(x) / length(x))
}
