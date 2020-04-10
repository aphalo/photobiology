#' Calculate spectral values at a different set of wavelengths
#'
#' Interpolate/re-express spectral irradiance (or other spectral quantity)
#' values at new wavelengths values. This is a low-level function operating
#' on numeric vectors and called by higher level functions in the package,
#' such as mathematical operators for classes for spectral data.
#'
#' @param w.length.in numeric vector of wavelengths (nm).
#' @param s.irrad a numeric vector of spectral values.
#' @param w.length.out numeric vector of wavelengths (nm).
#' @param fill a value to be assigned to out of range wavelengths.
#' @param ... additional arguments passed to \code{spline()}.
#'
#' @seealso \code{\link[stats]{splinefun}}.
#'
#' @return a numeric vector of interpolated spectral values.
#'
#' @export
#'
#' @note The current version of interpolate uses \code{spline} if fewer than 25
#' data points are available. Otherwise it uses \code{approx}. In the first case
#' a cubic spline is used, in the second case linear interpolation, which should
#' be faster.
#'
#' @family low-level functions operating on numeric vectors.
#'
#' @examples
#'
#' my.w.length <- 300:700
#' with(sun.data, interpolate_spectrum(w.length, s.e.irrad, my.w.length))
#'
interpolate_spectrum <- function(w.length.in,
                                 s.irrad,
                                 w.length.out,
                                 fill = NA,
                                 ...) {
  if (length(w.length.in) == length(w.length.out) &&
      all(w.length.in == w.length.out)) {
    # nothing to do
    return(s.irrad)
  }
  if (is.null(fill) && (w.length.out[1] < w.length.in[1] ||
                        w.length.out[length(w.length.out)] > w.length.in[length(w.length.in)])) {
    stop("Extrapolation attempted with fill == NULL")
  }
  selector <- w.length.out >= w.length.in[1] & w.length.out <= w.length.in[length(w.length.in)]
  s.irrad.out <- numeric(length(w.length.out))
  if (!is.null(fill)){
    s.irrad.out[!selector] <- fill
  }
  if (sum(selector) < 1) {
    NULL
  } else if (sum(selector) <= 25) {
    s.irrad.out[selector] <- stats::spline(x = w.length.in,
                                           y = s.irrad,
                                           xout = w.length.out[selector],
                                           ...)$y
  } else {
    s.irrad.out[selector] <- stats::approx(x = w.length.in,
                                           y = s.irrad,
                                           xout = w.length.out[selector])$y
  }
  return(s.irrad.out)
}
