#' Product from collection of spectra
#'
#' A method to compute the product of values across members of a collections of
#' spectra. Computes the product at each wavelength across all the spectra in the
#' collection returning a spectral object.
#'
#' @param x An R object. Currently this package defines methods for collections of
#'    spectral objects.
#' @param na.rm	logical. A value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param ...	Further arguments passed to or from other methods.
#'
#' @return If \code{x} is a collection spectral of objects, such as a
#'   "filter_mspct" object, the returned object is of same class as the
#'   members of the collection, such as "filter_spct", containing the product
#'   of the spectra.
#'
#' @note Omission of NAs is done separately at each wavelength. Interpolation is
#'   not applied, so all spectra in \code{x} must share the same set of
#'   wavelengths.
#'
#'   A product of spectral irradiance or spectral response is no longer a well
#'   defined physical quanttiy, and these product operations return an object of
#'   class generic_spct.
#'
#'   Objects of classes raw_spct and cps_spct can contain data from multiple
#'   scans. This functions are implemented for these classes only for the case
#'   when all member spectra contain data for a single scan, or spliced into a
#'   single column in the case of cps_spct members.
#'
#' @seealso See \code{\link[base]{prod}} for the \code{prod()} method used for
#'   the computations.
#'
#' @export
#'
s_prod <- function(x, na.rm, ...)
  UseMethod("s_prod")

#' @describeIn s_prod
#'
#' @export
#'
s_prod.default <- function(x, na.rm = FALSE, ...) {
  warning("Metod 's_prod()' not implementd for objects of class ",
          class(x)[1],
          ".")
  ifelse(is.any_mspct(x), do.call(class(x[[1]])[1], args = list()), NA)
}

#' @describeIn s_prod
#'
#' @export
#'
s_prod.source_mspct <- function(x, na.rm = FALSE, ...) {
  warning("A product of irradiance values does not yield response irradiance!!")
  rowwise_source(
    x,
    .fun = base::prod,
    na.rm = na.rm,
    col.name.tag = ".prod",
    .fun.name = "Product of"
  )
}

#' @describeIn s_prod
#'
#' @export
#'
s_prod.response_mspct <- function(x, na.rm = FALSE, ...) {
  warning("A product of response values does not yield response values!!")
  rowwise_response(
    x,
    .fun = base::prod,
    na.rm = na.rm,
    col.name.tag = ".prod",
    .fun.name = "Product of"
  )
}

#' @describeIn s_prod
#'
#' @export
#'
s_prod.filter_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_filter(x,
                 .fun = base::prod,
                 na.rm = na.rm,
                 .fun.name = "Product of")
}

#' @describeIn s_prod
#'
#' @export
#'
s_prod.reflector_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_reflector(x,
                    .fun = base::prod,
                    na.rm = na.rm,
                    .fun.name = "Product of")
}

#' @describeIn s_prod
#'
#' @export
#'
s_prod.calibration_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_calibration(
    x,
    .fun = base::prod,
    na.rm = na.rm,
    col.name.tag = ".prod",
    .fun.name = "Product of"
  )
}

#' @describeIn s_prod
#'
#' @export
#'
s_prod.cps_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_cps(
    x,
    .fun = base::prod,
    na.rm = na.rm,
    col.name.tag = ".prod",
    .fun.name = "Product of"
  )
}

#' @describeIn s_prod
#'
#' @export
#'
s_prod.raw_mspct <- function(x, na.rm = FALSE, ...) {
  rowwise_raw(
    x,
    .fun = base::prod,
    na.rm = na.rm,
    col.name.tag = ".prod",
    .fun.name = "Product of"
  )
}
