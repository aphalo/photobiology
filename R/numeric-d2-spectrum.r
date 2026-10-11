#' Calculate deuterium lamp output spectrum from fitted constants
#'
#' @description Calculate values by means of a nth degree polynomial from
#' user-supplied constants (for example from a lamp calibration certificate).
#'
#' @param w.length numeric vector of wavelengths (nm) for output
#' @param k a polynom:polynomial object with n constants for the polynomial
#' @param fill if NA, no extrapolation is done, and NA is returned for
#'   wavelengths outside the range 190 nm to 450 nm. If NULL then the tails are
#'   deleted. If 0 then the tails are set to zero, etc. NA is default.
#'
#' @return a dataframe with four numeric vectors with wavelength values
#'   (w.length), energy and photon irradiance (s.e.irrad, s.q.irrad) depending
#'   on the argument passed to unit.out (s.irrad).
#'
#' @export
#'
#' @note This is function is valid for wavelengths in the range 180 nm to 495
#'   nm, for wavelengths outside this range NAs are returned.
#' @examples
#' D2_spectrum(200)
#' D2_spectrum(170:220)
#'
D2_spectrum <- function(w.length,
                        k = photobiology::D2.UV653,
                        fill = NA_real_) {
  stopifnot(polynom::is.polynomial(k))
  poly.fun <- as.function(k)
  if (!is.null(fill)) {
    s.e.irrad <- ifelse(w.length >= 190 & w.length <= 450,
                        poly.fun(w.length),
                        fill)
  } else {
    w.length <- clip_wl(w.length, range = c(190, 450))
    s.e.irrad <- poly.fun(w.length)
  }
  out.data <- source_spct(w.length, s.e.irrad)
  comment(out.data) <- paste("Fitted spectrum for:", comment(k))
  return(out.data * 1e4)
}
