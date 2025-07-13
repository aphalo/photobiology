#' Mean from collection of spectra
#'
#' Method to compute the "parallel" mean of values across members of a
#' collection of spectra or of a spectral object containing multiple spectra in
#' long form.
#'
#' @details Method specializations compute the mean at each wavelength across a
#'   group of spectra stored in an object of one of the classes defined in
#'   package 'photobiolgy'. Trimming of extreme values is
#'   possible (trimmed mean) and omission of NAs is done separately at each
#'   wavelength. Interpolation is not applied, so all spectra in \code{x} must
#'   share the same set of wavelengths. An error is triggered if this condition
#'   is nor fulfilled.
#'
#' @param x An R object.
#' @param trim	numeric The fraction (0 to 0.5) of observations to be trimmed
#'   from each end of x before the mean is computed. Values of trim outside
#'   that range are taken as the nearest endpoint.
#' @param na.rm	logical A value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param ...	Further arguments passed to or from other methods.
#'
#' @return If \code{x} is a collection spectral of objects, such as a
#'   \code{"filter_mspct"} object, the returned object belongs to the same class
#'   as the members of the collection, such as \code{"filter_spct"}, containing
#'   the summary spectrum, with variables with names tagged for summaries other
#'   than mean or median.
#'
#' @note Objects of classes \code{raw_spct} and \code{cps_spct} can contain data
#'   from multiple scans in multiple variables or "columns". The parallel
#'   summaries' methods accept as arguments objects of these classes only if
#'   spectra contain data for a single spectrometer scan. In the case of
#'   \code{cps_spct} objects, a single column can also contain data from
#'   multiple scans spliced into a single variable.
#'
#' @section Deepest Curves: Parallel summaries differ fundamentally from the
#'   "deepest curves" obtained through functional data analysis (FDA) in that in
#'   functional data analysis one of the input curves is returned as the deepest
#'   one based on a decision criterion. In contrast the parallel summaries from
#'   package 'photobioloy' return one or more "fictional" curves different to
#'   any of those passed as inputs. This curve is constructed from independent
#'   summaries at each wavelength value.
#'
#' @seealso See \code{\link[base]{mean}} for the \code{mean()} method used for
#'   the computations.
#'
#' @export
#'
#' @examples
#' s_mean(sun_evening.mspct)
#'
s_mean <- function(x, trim, na.rm, ...) UseMethod("s_mean")

#' @rdname s_mean
#'
#' @export
#'
s_mean.default <- function(x, trim = 0, na.rm = FALSE, ...) {
  warning("Metod 's_mean()' not implementd for objects of class ",
          class(x)[1], ".")
  ifelse(is.any_mspct(x), do.call(class(x[[1]])[1], args = list()), NA)
}

#' @rdname s_mean
#'
#' @export
#'
s_mean.generic_spct <- function(x, trim = 0, na.rm = FALSE, ...) {
  if (getMultipleWl(x) > 1) {
    s_mean(subset2mspct(x), trim = trim, na.rm = na.rm, ...)
  } else {
    x
  }
}

#' @rdname s_mean
#'
#' @export
#'
s_mean.source_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {
  rowwise_source(x,
                 .fun = base::mean,
                 trim = trim,
                 na.rm = na.rm,
                 .fun.name = "Mean of")
}

#' @rdname s_mean
#'
#' @export
#'
s_mean.response_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {
  rowwise_response(x,
                   .fun = base::mean,
                   trim = trim,
                   na.rm = na.rm,
                   .fun.name = "Mean of")
}

#' @rdname s_mean
#'
#' @export
#'
s_mean.filter_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {
  rowwise_filter(x,
                 .fun = base::mean,
                 trim = trim,
                 na.rm = na.rm,
                 .fun.name = "Mean of")
}

#' @rdname s_mean
#'
#' @export
#'
s_mean.reflector_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {
  rowwise_reflector(x,
                    .fun = base::mean,
                    trim = trim,
                    na.rm = na.rm,
                    .fun.name = "Mean of")
}

#' @rdname s_mean
#'
#' @export
#'
s_mean.calibration_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {
  rowwise_calibration(x,
                      .fun = base::mean,
                      trim = trim,
                      na.rm = na.rm,
                      .fun.name = "Mean of")
}

#' @rdname s_mean
#'
#' @export
#'
s_mean.cps_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {
  rowwise_cps(x,
              .fun = base::mean,
              trim = trim,
              na.rm = na.rm,
              .fun.name = "Mean of")
}

#' @rdname s_mean
#'
#' @export
#'
s_mean.raw_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {
  rowwise_raw(x,
              .fun = base::mean,
              trim = trim,
              na.rm = na.rm,
              .fun.name = "Mean of")
}
