#' Mean plus and minus standard error from collection of spectra
#'
#' Method to compute the "parallel" mean and limits based on SEM. The spectral
#' values are summarised across members of a collection of spectra or of a
#' spectral object containing multiple spectra in long form.
#'
#' @details Method specializations compute the mean and limits based on SEM at
#'   each wavelength across a group of spectra stored in an object of one of the
#'   classes defined in package 'photobiology'. Omission of NAs is done
#'   separately at each wavelength. Interpolation is not applied, so all spectra
#'   in \code{x} must share the same set of wavelengths. An error is triggered
#'   if this condition is nor fulfilled. The value passed as argument to `mult`
#'   can be used to estimate a confidence interval for each mean value.
#'
#' @inheritParams s_mean
#' @param mult numeric number of multiples of standard error.
#'
#' @inherit s_mean note return
#'
#' @inheritSection s_mean Deepest Curves
#'
#' @seealso See \code{\link[base]{mean}} for the \code{mean()} method used for
#'   the computations.
#'
#' @export
#'
#' @examples
#' s_mean_se_band(sun_evening.mspct)
#'
s_mean_se_band <- function(x, na.rm, mult, ...)
  UseMethod("s_mean_se_band")

#' @rdname s_mean_se_band
#'
#' @export
#'
s_mean_se_band.default <- function(x, na.rm = FALSE, mult = 1, ...) {
  warning("Metod 's_mean_se_band()' not implementd for objects of class ",
          class(x)[1],
          ".")
  ifelse(is.any_mspct(x), do.call(class(x[[1]])[1], args = list()), NA)
}

#' @rdname s_mean_se_band
#'
#' @export
#'
s_mean_se_band.generic_spct <- function(x, na.rm = FALSE, mult = 1, ...) {
  s_mean_se_band(subset2mspct(x), na.rm = na.rm, mult = mult, ...)
}

#' @rdname s_mean_se_band
#'
#' @export
#'
s_mean_se_band.filter_mspct <-
  function(x, na.rm = FALSE, mult = 1, ...) {
    rowwise_filter(
      x,
      .fun = list(base::mean, se.m, se.p),
      col.name.tag = c("", ".se.m", ".se.p"),
      na.rm = na.rm,
      mult = mult,
      .fun.name = "Mean and mean +/- SEM of"
    )
  }

#' @rdname s_mean_se_band
#'
#' @export
#'
s_mean_se_band.source_mspct <-
  function(x, na.rm = FALSE, mult = 1, ...) {
    rowwise_source(
      x,
      .fun = list(base::mean, se.m, se.p),
      col.name.tag = c("", ".se.m", ".se.p"),
      na.rm = na.rm,
      mult = mult,
      .fun.name = "Mean and mean +/- SEM of"
    )
  }

#' @rdname s_mean_se_band
#'
#' @export
#'
s_mean_se_band.response_mspct <-
  function(x, na.rm = FALSE, mult = 1, ...) {
    rowwise_response(
      x,
      .fun = list(base::mean, se.m, se.p),
      col.name.tag = c("", ".se.m", ".se.p"),
      na.rm = na.rm,
      mult = mult,
      .fun.name = "Mean and mean +/- SEM of"
    )
  }

#' @rdname s_mean_se_band
#'
#' @export
#'
s_mean_se_band.reflector_mspct <-
  function(x, na.rm = FALSE, mult = 1, ...) {
    rowwise_reflector(
      x,
      .fun = list(base::mean, se.m, se.p),
      col.name.tag = c("", ".se.m", ".se.p"),
      na.rm = na.rm,
      mult = mult,
      .fun.name = "Mean and mean +/- SEM of"
    )
  }

#' @rdname s_mean_se_band
#'
#' @export
#'
s_mean_se_band.calibration_mspct <-
  function(x, na.rm = FALSE, mult = 1, ...) {
    rowwise_calibration(
      x,
      .fun = list(base::mean, se.m, se.p),
      col.name.tag = c("", ".se.m", ".se.p"),
      na.rm = na.rm,
      mult = mult,
      .fun.name = "Mean and mean +/- SEM of"
    )
  }

#' @rdname s_mean_se_band
#'
#' @export
#'
s_mean_se_band.cps_mspct <- function(x, na.rm = FALSE, mult = 1, ...) {
  rowwise_cps(
    x,
    .fun = list(base::mean, se.m, se.p),
    col.name.tag = c("", ".se.m", ".se.p"),
    na.rm = na.rm,
    mult = mult,
    .fun.name = "Mean and mean +/- SEM of"
  )
}

#' @rdname s_mean_se_band
#'
#' @export
#'
s_mean_se_band.raw_mspct <- function(x, na.rm = FALSE, mult = 1, ...) {
  rowwise_raw(
    x,
    .fun = list(base::mean, se.m, se.p),
    col.name.tag = c("", ".se.m", ".se.p"),
    na.rm = na.rm,
    mult = mult,
    .fun.name = "Mean and mean +/- SEM of"
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
se.p <- function(x, na.rm = FALSE, mult = 1, ...) {
  stopifnot(is.numeric(x))
  if (na.rm) {
    x <- stats::na.omit(x)
  }
  mean(x) + mult * sqrt(stats::var(x) / length(x))
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
se.m <- function(x, na.rm = FALSE, mult = 1, ...) {
  stopifnot(is.numeric(x))
  if (na.rm) {
    x <- stats::na.omit(x)
  }
  mean(x) - mult * sqrt(stats::var(x) / length(x))
}
