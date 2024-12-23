#' Extraterrestrial irradiance
#'
#' Estimate of down-welling solar (short wave) irradiance at the top of the
#' atmosphere above a location on Earth, computed based on angles, Sun-Earth
#' distance and the solar constant. Astronomical computations are done with
#' function \code{sun_angles()}.
#'
#' @param time A "vector" of POSIXct Time, with any valid time zone (TZ) is
#'   allowed, default is current time.
#' @param tz character string indicating time zone to be used in output.
#' @param geocode data frame with variables lon and lat as numeric values
#'   (degrees), nrow > 1, allowed.
#' @param solar.constant numeric or character If character, "WMO" or "NASA", if
#'   numeric, an irradiance value in the same units as the value to be returned.
#'
#' @return Numeric vector of extraterrestrial irradiance (in W / m2 if solar
#'   constant is a character value).
#'
#' @seealso Function \code{\link{sun_angles}}.
#'
#' @examples
#' library(lubridate)
#'
#' irrad_extraterrestrial(ymd_hm("2021-06-21 12:00", tz = "UTC"))
#'
#' irrad_extraterrestrial(ymd_hm("2021-12-21 20:00", tz = "UTC"))
#'
#' irrad_extraterrestrial(ymd_hm("2021-06-21 00:00", tz = "UTC") + hours(1:23))
#'
#' @export
#'
irrad_extraterrestrial <-
  function(time = lubridate::now(tzone = "UTC"),
           tz = lubridate::tz(time),
           geocode = tibble::tibble(lon = 0, lat = 51.5, address = "Greenwich"),
           solar.constant = "NASA") {
    solar.cnst.map <- c(NASA = 1360, WMO = 1367) # W / m2
    if (is.character(solar.constant)) {
      solar.constant <- solar.cnst.map[solar.constant]
    }
    angles <- sun_angles(time = time,
                         tz = tz,
                         geocode = geocode,
                         use.refraction = FALSE) # no atmosphere!
    rel.distance <- angles[["distance"]]
    sun.elevation <- angles[["elevation"]]
    ifelse(sun.elevation <= 0,
           0,
           solar.constant *
             cos((90 - sun.elevation) / 180 * pi) / # degrees -> radians
             rel.distance^2)
  }
