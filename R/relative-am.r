#' Relative Air Mass (AM)
#'
#' Approximate relative air mass (AM) from sun elevation or
#' sun zenith angle.
#'
#' @param elevation.angle,zenith.angle numeric vector Angle in degrees for the
#'   sun position. An argument should be passed to one and only one of
#'   \code{elevation_angle} or \code{zenith_angle}.
#' @param occluded.value numeric Value to return when elevation angle is
#'   negative (sun below the horizon).
#'
#' @details This is an implementation of equation (3) in Kasten and Young
#'   (1989). This equation is only an approximation to the tabulated values in
#'   the same paper. Returned values are rounded to three significant digits.
#'
#' @note Although relative air mass is not defined when the sun is not visible,
#'   returning a value different from the default \code{NA} might be useful in
#'   some cases.
#'
#' @export
#'
#' @references
#' F. Kasten, A. T. Young (1989) Revised optical air mass tables and
#' approximation formula. Applied Optics, 28, 4735-. doi:10.1364/ao.28.004735.
#'
#' @examples
#'
#' relative_AM(c(90, 60, 30, 1, -10))
#' relative_AM(c(90, 60, 30, 1, -10), occluded.value = Inf)
#' relative_AM(zenith.angle = 0)
#'
relative_AM <- function(elevation.angle = NULL,
                        zenith.angle = NULL,
                        occluded.value = NA) {
  stopifnot(xor(is.null(elevation.angle), is.null(zenith.angle)))
  if (is.null(elevation.angle)) {
    elevation.angle <- 90 - zenith.angle
  }
  stopifnot(all(elevation.angle >= -90 & elevation.angle <= 90))
  signif(
    ifelse(elevation.angle > 0,
           (sin(elevation.angle * pi / 180) + (0.1500 * (elevation.angle + 3.885)^-1.253))^-1,
           occluded.value[1]),
    3)
}
