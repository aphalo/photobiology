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
#' @param fill a numeric value to be assigned to out of range wavelengths.
#' @param method character string One of \code{"auto"}, \code{"approx"},
#'   \code{"spline"}, \code{"skip"}.
#' @param ... additional arguments passed to \code{spline()}.
#'
#' @seealso \code{\link[stats:splinefun]{spline}()} and
#'   \code{\link[stats:approxfun]{approx}()}.
#'
#' @return a numeric vector of interpolated spectral values.
#'
#' @export
#'
#' @details Depending on \code{method} natural spline interpolation or linear
#'   interpolation are used. With \code{method = spline} a call to
#'   \code{\link[stats:splinefun]{spline}} with \code{method = "natural"} is
#'   used and with \code{method = "approx"} a call to
#'   \code{\link[stats:approxfun]{approx}} is used. If \code{method = "auto"} or
#'   \code{method = NULL} when 100 or fewer distinct wavelengths are available
#'   as input and/or the maximum wavelength step size in \code{w.length.in} is
#'   more than three times the minimum wavelength step size in
#'   \code{w.length.out} \code{"spline"} is used and \code{"approx"} otherwise.
#'   Finally, with \code{method = "skip"} the input is returned unchanged.
#'
#'   If \code{w.length.out} is a numeric vector and \code{length.out = NULL}, it
#'   directly gives the target wavelengths for interpolation. If it is
#'   \code{NULL}, and \code{length.out} is an integer value evenly spaced
#'   wavelength values covering the same wavelength range as in the input are
#'   generated. If \code{w.length.out} is a numeric vector and \code{length.out}
#'   is an integer value, \code{length.out} evenly spaced wavelengths covering
#'   the wavelength range of \code{w.length.out} are generated.
#'   \emph{Extrapolation is not supported.}
#'
#'   With default \code{fill = NA} if the output exceeds the wavelength range of
#'   the input, extrapolated values are filled with \code{NA} values. With
#'   \code{fill = NULL} wavelengths outside the wavelength range of input data
#'   are discarded. A numerical value can be also be provided as fill. While
#'   \code{interpolate_spectrum} supports interpolation of a single numeric
#'   vector, \code{interpolate_wl} applies, one at a time, interpolation to all
#'   numeric columns found in \code{x}.
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
                                 method = "approx",
                                 ...) {
  if (length(w.length.in) == length(w.length.out) &&
        all(w.length.in == w.length.out)) {
    # nothing to do
    return(s.irrad)
  }
  if (is.null(fill) &&
        (w.length.out[1] < w.length.in[1] ||
           w.length.out[length(w.length.out)] >
             w.length.in[length(w.length.in)])) {
    stop("Extrapolation attempted with fill == NULL")
  }
  selector <- w.length.out >= w.length.in[1] &
                w.length.out <= w.length.in[length(w.length.in)]
  if (sum(selector) < 1 || length(w.length.in) < 2) {
    method <- "skip"
  }
  step.size.ratio <-
    photobiology::stepsize(w.length.out)[1] /
    photobiology::stepsize(w.length.in)[2]
  if (is.null(method) || method == "auto") {
    if (length(w.length.in) <= 100 || step.size.ratio > 3) {
      method <- "spline"
    } else {
      method <- "approx"
    }
  }
  s.irrad.out <- numeric(length(w.length.out))
  if (!is.null(fill)) {
    s.irrad.out[!selector] <- fill
  }
  if (method == "spline") {
    s.irrad.out[selector] <-
      stats::spline(x = w.length.in,
                    y = s.irrad,
                    xout = w.length.out[selector],
                    method = "natural",
                    ...)[["y"]]
  } else if (method == "approx") {
    s.irrad.out[selector] <-
      stats::approx(x = w.length.in,
                    y = s.irrad,
                    xout = w.length.out[selector])[["y"]]
  } else if (method != "skip") {
    stop("Wrong 'method' argument: ", method)
  }
  return(s.irrad.out)
}
