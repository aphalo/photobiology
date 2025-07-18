#' Sum from collection of spectra
#'
#' Method to compute the "parallel" sum of values across members of a
#' collection of spectra or of a spectral object containing multiple spectra in
#' long form.
#'
#' @details Method specializations compute the sum at each wavelength across a
#'   group of spectra stored in an object of one of the classes defined in
#'   package 'photobiology'. Omission of NAs is done
#'   separately at each wavelength. Interpolation is not applied, so all spectra
#'   in \code{x} must share the same set of wavelengths. An error is triggered
#'   if this condition is nor fulfilled.
#'
#' @inheritParams s_mean
#'
#' @inherit s_mean note return
#'
#' @note The sum operation is meaningful only for certain physical
#'   quantities or bases of expression.
#'
#' @seealso See \code{\link[base]{sum}} for the \code{sum()} method used for
#'   the computations.
#'
#' @export
#'
#' @examples
#' s_sum(sun_evening.mspct)
#'
s_sum <- function(x, na.rm, ...)
  UseMethod("s_sum")

#' @rdname s_sum
#'
#' @export
#'
s_sum.default <- function(x, na.rm = FALSE, ...) {
  warning("Metod 's_sum()' not implementd for objects of class ",
          class(x)[1],
          ".")
  ifelse(is.any_mspct(x), do.call(class(x[[1]])[1], args = list()), NA)
}

#' @rdname s_sum
#'
#' @export
#'
s_sum.generic_spct <- function(x, na.rm = FALSE, ...) {
  s_sum(subset2mspct(x), na.rm = na.rm, ...)
}

#' @rdname s_sum
#'
#' @export
#'
s_sum.filter_mspct <- function(x, na.rm = FALSE, ...) {
  warning("A sum of Tfr values does not yield Tfr values, ",
          "while a summ of A values yields A values!!")
  rowwise_filter(
    x,
    .fun = base::sum,
    na.rm = na.rm,
    col.name.tag = ".sum",
    .fun.name = "Sum of"
  )
}

#' @rdname s_sum
#'
#' @export
#'
s_sum.source_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_source(x,
                 .fun = base::sum,
                 na.rm = na.rm,
                 .fun.name = "Sum of")
}

#' @rdname s_sum
#'
#' @export
#'
s_sum.response_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_response(x,
                   .fun = base::sum,
                   na.rm = na.rm,
                   .fun.name = "Sum of")
}

#' @rdname s_sum
#'
#' @export
#'
s_sum.reflector_mspct <- function(x, na.rm = FALSE, ...) {
  warning("A sum of Rfr values does not yield Rfr values!!")
  rowwise_reflector(
    x,
    .fun = base::sum,
    na.rm = na.rm,
    col.name.tag = ".sum",
    .fun.name = "Sum of"
  )
}

#' @rdname s_sum
#'
#' @export
#'
s_sum.calibration_mspct <- function(x, na.rm = FALSE, ...) {
  warning("A sum of irrad.mult values does not yield irrad.mult values!!")
  rowwise_calibration(
    x,
    .fun = base::sum,
    na.rm = na.rm,
    col.name.tag = ".sum",
    .fun.name = "Sum of"
  )
}

#' @rdname s_sum
#'
#' @export
#'
s_sum.cps_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_cps(
    x,
    .fun = base::sum,
    na.rm = na.rm,
    col.name.tag = ".sum",
    .fun.name = "Sum of"
  )
}

#' @rdname s_sum
#'
#' @export
#'
s_sum.raw_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_raw(
    x,
    .fun = base::sum,
    na.rm = na.rm,
    col.name.tag = ".sum",
    .fun.name = "Sum of"
  )
}
