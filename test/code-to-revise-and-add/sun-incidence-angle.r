#' Angle of incidence of the direct solar beam
#'
#' The angle of incidence of the direct solar beam on a plane oriented to an
#' arbitrary azimuth and with an arbitrary tilt with respect to the horizontal
#' are computed either from solar angles passed as argument or by first
#' computing these angles.
#'
#' @param time A "vector" of POSIXct Time, with any valid time zone (TZ) is
#'   allowed, default is current time.
#' @param tz character string indicating time zone to be used in output.
#' @param geocode data frame with variables lon and lat as numeric values
#'   (degrees), nrow > 1, allowed.
#' @param use.refraction logical Flag indicating whether to correct for
#'   fraction in the atmosphere.
#' @param plane.azimuth The azimuth angle in degrees of the plane surface
#'   receiving the direct solar beam, measured in the same way as the
#'   solar azimuth..
#' @param plane.tilt The tilt angle in degrees between the normal to the
#'   plane surface receiving the direct solar beam and the normal to a
#'   horizontal plane.
#' @param sun.angles The output of function \code{sun_angles()} or a data.frame
#'   with at least columns named \code{azimuth} and \code{elevation} containing
#'   values expressed in degrees. If \code{NULL} these are calculated by calling
#'   \code{sun_angles()}.
#'
#' @return Numeric vector of angles in degrees, with 90 degrees indicating
#'   normal incidence of the beam. The sign is removed and
#'   if the beam hits the back of the plane, \code{Inf} is returned.
#'
#' @examples
#'
#' sun_incidence_angle(plane.azimuth = 90, # East facing butirrelevant
#'                     plane.tilt = 0, # horizontal
#'                     sun.angles = data.frame(azimuth = 90,
#'                                             elevation = 45))
#'
#' sun_incidence_angle(plane.azimuth = 90, # East facing
#'                     plane.tilt = 45,
#'                     sun.angles = data.frame(azimuth = 90,
#'                                             elevation = 45))
#'
#' sun_incidence_angle(plane.azimuth = 90, # West facing
#'                     plane.tilt = 45,
#'                     sun.angles = data.frame(azimuth = 0,
#'                                             elevation = 0))
#' @export
#'
sun_incidence_angle <- function(time = lubridate::now(),
                                tz = lubridate::tz(time),
                                geocode = tibble::tibble(lon = 0,
                                                         lat = 51.5,
                                                         address = "Greenwich"),
                                use.refraction = FALSE,
                                plane.azimuth = 0,
                                plane.tilt = 0,
                                plane.back.value = NA_real_,
                                sun.angles = NULL) {
  if (is.null(sun.angles)) {
    sun.angles <-
      sun_angles(time = time,
                 tz = tz,
                 geocode = geocode,
                 use.refraction = use.refraction)

  }
  sun.zenith.angle <- (90 - sun.angles[["elevation"]]) * pi / 180
  # is the sun below the horizon
  sun.zenith.angle <- ifelse(sun.zenith.angle > pi / 2, NA_real_, sun.zenith.angle)
  sun.azimuth <- (sun.angles[["azimuth"]]) * pi / 180
  plane.tilt <- plane.tilt * pi / 180
  plane.azimuth <- plane.azimuth * pi / 180
  # the sign of the angle is lost in this equation
  incidence <-
    acos(cos(sun.zenith.angle) * cos(plane.tilt) +
           sin(plane.tilt) * sin(sun.zenith.angle) *
           cos(sun.azimuth + plane.azimuth)) * 180 / pi
  delta.azimuth <- ifelse(sun.azimuth >= plane.azimuth,
                          sun.azimuth - plane.azimuth,
                          plane.azimuth - sun.azimuth)
  # check sign, is the plane illuminated from the back?
  incidence <- ifelse(delta.azimuth > pi / 2,
                      -incidence,
                      incidence)
  ifelse(is.null(plane.back.value) & incidence < 0,
         plane.back.value,
         incidence)
}

