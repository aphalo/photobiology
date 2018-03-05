#' Relative Air Mass (AM)
#'
#' Approximate relative air mass (AM) from sun elevation or
#' sun zenith angle.
#'
#' @param elevation_angle,zenith_angle numeric vector Angle in degrees for the
#' sun position. An argument should be passed to one and only one of \code{elevation_angle} or
#' \code{zenith_angle}.
#'
#' @detail This is an implementation of equation (3) in Kasten and Young (1989).
#'
#' @export
#'
#' @references
#' F. Kasten, A. T. Young (1989) Revised optical air mass tables and
#' approximation formula. Applied Optics, 28, 4735-. doi:10.1364/ao.28.004735.
#'
relative_AM <- function(elevation_angle = NULL, zenith_angle = NULL) {
  stopifnot(xor(is.null(elevation), is.null(zenith_angle)))
  elevation_angle <- ifelse(is.null(elevation_angle), 90 - zenith_angle, elevation_angle)
  stopifnot(elevation_angle >= -90 || elevation_angle <= 90)
  ifelse(elevation_angle >= 0,
         (sin(elevation_angle * pi / 180) + (0.1500 * (elevation_angle + 3.885)^-1.253))^-1,
         NA_real_)
}
