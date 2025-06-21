#' Reflectance at a planar boundary
#'
#' The reflectance at the planar boundary between two media, or interface, can
#' be computed from the relative refractive index. Reflectance depends on
#' polarization, and the process of reflection can generate polarized light
#' through selective reflection of \eqn{s} and \eqn{p} components. A perfectly
#' flat (i.e., polished) interface creates specular reflection, and this is the
#' case that these functions describe. These function describe a single
#' interface, and for example in a glass pane, a light beam will cross two
#' air-glass interfaces.
#'
#' @param angle_deg,angle numeric vector Angle of incidence of the light beam,
#'   in degrees or radians. If both are supplied, radians take precedence.
#' @param n numeric vector, or generic_spct object Relative refractive index.
#'   The default 1.5 is suitable for crown glass or acrylic interacting with
#'   visible light. \eqn{n} depends on wavelength, more or less strongly
#'   depending on the material.
#' @param p_fraction numeric in range 0 to 1. Polarization, defaults to 0.5
#'   assuming light that is not polarized.
#'
#' @details These functions implement Fresnel's formulae. All parameters accept
#'   vectors as arguments. If both n and angle are vectors with length different
#'   from one, they should both have the same length. Reflectance depends on
#'   polarization, the \eqn{s} and \eqn{p} components need to be computed
#'   separately and added up. \code{Rfr_from_n()} is for non-polarized light,
#'   i.e., with equal contribution of the two components.
#'
#' @return If \code{n} is a numeric vector the returned value is a vector of
#'   reflectances, while if \code{n} is a \code{generic_spct} object the
#'   returned value is a \code{reflector_spct} object.
#'
#' @export
#'
#' @examples
#'
#' Rfr_from_n(0:90)
#' Rfr_from_n(0:90, p_fraction = 1)
#' Rfr_from_n(0:90, n = 1.333) # water
#'
Rfr_from_n <- function(angle_deg,
                       angle = angle_deg / 180 * pi,
                       n = 1.5,
                       p_fraction = 0.5) {
  stopifnot(all(p_fraction >= 0 & p_fraction <= 1))
  stopifnot(all(angle >= 0 & angle <= pi / 2))
  if (is.generic_spct(n)) {
    z <- reflector_spct(w.length = n[["w.length"]],
                        Rfr = Rfr_from_n(angle = angle,
                                         n = n[["n"]], p_fraction = p_fraction),
                        comment = paste("Computed from:", what_measured(n)),
                        Rfr.type = "total")
  } else {
    z <-  Rfr_p_from_n(angle = angle, n = n) * p_fraction +
      Rfr_s_from_n(angle = angle, n = n) * (1 - p_fraction)
  }
  z
}

#' @rdname Rfr_from_n
#'
#' @export
#'
Rfr_p_from_n <- function(angle_deg, angle = angle_deg / 180 * pi, n = 1.5) {
  (
    (n^2 * cos(angle) - sqrt(n^2  - sin(angle)^2)) /
      (n^2 * cos(angle) + sqrt(n^2  - sin(angle)^2))
  )^2
}

#' @rdname Rfr_from_n
#'
#' @export
#'
Rfr_s_from_n <- function(angle_deg, angle = angle_deg / 180 * pi, n = 1.5) {
  (
    (cos(angle) - sqrt(n^2  - sin(angle)^2)) /
      (cos(angle) + sqrt(n^2  - sin(angle)^2))
  )^2
}
