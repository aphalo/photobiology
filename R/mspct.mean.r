#' Mean from collection of spectra
#'
#' A method to compute the mean of values across members of a collections of
#' spectra. Computes the mean at each wavelength across all the spectra in the
#' collection returning a spectral object.
#'
#' @param x An R object. Currently this package defines methods for collections of
#'    spectral objects.
#' @param trim	numeric. The fraction (0 to 0.5) of observations to be trimmed from
#'   each end of x before the mean is computed. Values of trim outside that
#'   range are taken as the nearest endpoint.
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
#' each wavelength. Interpolation is not applied, so all spectra in \code{x}
#' must share the same set of wavelengths.
#'
#' @seealso See \code{\link[base]{mean}} for the \code{mean()} method used for
#'   the computations.
#'
#' @export
#'
s_mean <- function(x, trim, na.rm, ...) UseMethod("s_mean")

#' @describeIn s_mean
#'
#' @export
#'
s_mean.default <- function(x, trim = 0, na.rm = FALSE, ...) {
  warning("Metod 's_mean()' not implementd for objects of class ", class(x)[1], ".")
  ifelse(is.any_mspct(x), do.call(class(x[[1]])[1], args = list()), NA)
}

#' @describeIn s_mean
#'
#' @export
#'
s_mean.filter_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {
  rowwise_filter(x, .fun = base::mean, trim = trim, na.rm = na.rm, .fun.name = "Mean of")
}

#' @describeIn s_mean
#'
#' @export
#'
s_mean.source_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {
  rowwise_source(x, .fun = base::mean, trim = trim, na.rm = na.rm, .fun.name = "Mean of")
}

#' @describeIn s_mean
#'
#' @export
#'
s_mean.response_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {
  rowwise_response(x, .fun = base::mean, trim = trim, na.rm = na.rm, .fun.name = "Mean of")
}

#' @describeIn s_mean
#'
#' @export
#'
s_mean.reflector_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {
  rowwise_reflector(x, .fun = base::mean, trim = trim, na.rm = na.rm, .fun.name = "Mean of")
}
