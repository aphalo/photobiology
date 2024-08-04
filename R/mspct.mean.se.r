#' Mean and standard error from collection of spectra
#'
#' Method to compute the "parallel" mean and the SEM. The spectral
#' values are summarised across members of a collection of spectra or of a
#' spectral object containing multiple spectra in long form.
#'
#' @details Method specializations compute the  mean and SEM at
#'   each wavelength across a group of spectra stored in an object of one of the
#'   classes defined in package 'photobiology'. Omission of NAs is done
#'   separately at each wavelength. Interpolation is not applied, so all spectra
#'   in \code{x} must share the same set of wavelengths. An error is triggered
#'   if this condition is nor fulfilled. The value passed as argument to `mult`
#'
#' @param x An R object Currently this package defines methods for collections of
#'    spectral objects.
#' @param na.rm	logical A value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param mult	numeric number of multiples of standard error.
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
#' @seealso See \code{\link[base]{mean}} for the \code{mean()} method to
#'   compute the mean and \code{\link{se}} for the method used to compute the
#'   standard error of the mean.
#'
#' @export
#'
#' @examples
#' s_mean_se(sun_evening.mspct)
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
s_mean_se.generic_spct <- function(x, na.rm = FALSE, mult = 1, ...) {
  s_mean_se(subset2mspct(x), na.rm = na.rm, mult = mult, ...)
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
