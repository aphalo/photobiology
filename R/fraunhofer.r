#' Diffraction
#'
#' Diffraction of optical radiation passing through a single slit can
#' be computed with function \code{diffraction_single_slit()}, which
#' implements Fraunhofer's equation. Diffraction plus interference for a
#' pair of slits can be computed with \code{diffraction_double_slit()}.
#'
#' @param w.length numeric Wavelength (nm).
#' @param slit.width numeric Width of the slit (m).
#' @param angle numeric vector Angle (radians).
#'
#' @return A numeric vector of the same length as \code{angle}, containing
#'   relative intensities.
#'
#' @export
#'
#' @examples
#' diffraction_single_slit(w.length = 550,
#'                              slit.width = 1e-5,
#'                              angle = 0)
#'
#' # use odd number for length.out so that 0 is in the sequence
#' angles <- pi * seq(from = -1/2, to = 1/2, length.out = 501)
#'
#' plot(angles,
#'      diffraction_single_slit(w.length = 550, # 550 nm
#'                              slit.width = 6e-6, # 6 um
#'                              angle = angles),
#'      type = "l",
#'      ylab = "Relative irradiance (/1)",
#'      xlab = "Angle (radian)")
#'
#' plot(angles,
#'      diffraction_double_slit(w.length = 550, # 550 nm
#'                              slit.width = 6e-6, # 6 um
#'                              slit.distance = 18e-6, # 18 um
#'                              angle = angles),
#'      type = "l",
#'      ylab = "Relative irradiance (/1)",
#'      xlab = "Angle (radian)")
#'
diffraction_single_slit <- function(w.length,
                                    slit.width,
                                    angle) {
  w.length <- w.length * 1e-9 # nm -> m

  sinc(slit.width / w.length * sin(angle))^2
}

#' @rdname diffraction_single_slit
#'
#' @param slit.distance numeric Distance between the centres of the two slits
#'   (m).
#'
#' @export
#'
diffraction_double_slit <- function(w.length,
                                    slit.width,
                                    slit.distance,
                                    angle) {
  w.length <- w.length * 1e-9 # nm -> m
  a <- slit.width / 2
  d <- slit.distance / 2
  sinc(angle * a / w.length)^2 * cos(angle * d / w.length)^2
}

#' sinc
#'
#' @keywords internal
#'
sinc <- function(x) {
  ifelse(abs(x) < 1e-10, 1, sin(x) / x)
}



