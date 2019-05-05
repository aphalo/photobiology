#' Solar angles
#'
#' This function returns the solar angles at a given time and location.
#'
#' @param time A "vector" of POSIXct Time, with any valid time zone (TZ) is
#'   allowed, default is current time.
#' @param tz character string indicating time zone to be used in output.
#' @param geocode data frame with variables lon and lat as numeric values
#'   (degrees), nrow > 1, allowed.
#' @param use.refraction logical Flag indicating whether to correct for
#'   fraction in the atmosphere.
#'
#' @return A data.frame with variables time (in same TZ as input), TZ, solartime,
#'   longitude, latitude, address, azimuth, and elevation. If a data frame with
#'   multiple rows is passed to \code{geocode} and a vector of times longer
#'   than one is passed to \code{time}, sun position for all combinations of
#'   locations and times are returned are returned by \code{sun_angles}. In
#'   contrast, convenience functions returning a vector.
#'
#' @family astronomy related functions
#'
#' @details This function is an implementation of Meeus equations as used in
#'   NOAAs on-line web calculator, which are precise and valid for a very broad
#'   range of dates (years -1000 to 3000 at least). The apparent solar
#'   elevations near sunrise and sunset are affected by refraction in the
#'   atmosphere, which does in turn depend on weather conditions. The effect of
#'   refraction on the apparent position of the sun is only an estimate based on
#'   "typical" conditions for the atmosphere. The computation is not defined for
#'   latitudes 90 and -90 degrees, i.e. exactly at the poles.
#'
#'   In the current implementation functions \code{sun_azimuth},
#'   \code{sun_elevation}, and \code{sun_zenith_angle} are wrappers
#'   on \code{sun_angles}, so if more than one angle is needed it is
#'   preferable to directly call \code{sun_angles} as it will be faster.
#'
#' @note
#'   There exists a different R implementation of the same algorithms called
#'   "AstroCalcPureR" available as function \code{astrocalc4r} in package
#'   'fishmethods'. Although the equations used are almost all the same, the
#'   function signatures and which values are returned differ. In particular,
#'   the present implementation splits the calculation into two separate
#'   functions, one returning angles at given instants in time, and a
#'   separate one returning the timing of events for given dates.
#'
#' @export
#' @examples
#' library(lubridate)
#' sun_angles()
#' sun_azimuth()
#' sun_elevation()
#' sun_zenith_angle()
#' sun_angles(ymd_hms("2014-09-23 12:00:00"))
#' sun_angles(ymd_hms("2014-09-23 12:00:00"),
#'            geocode = data.frame(lat=60, lon=0))
#' sun_angles(ymd_hms("2014-09-23 12:00:00") + minutes((0:6) * 10))
#'
sun_angles <- function(time = lubridate::now(tzone = "UTC"),
                       tz = lubridate::tz(time),
                       geocode = tibble::tibble(lon = 0,
                                                lat = 51.5,
                                                address = "Greenwich"),
                       use.refraction = FALSE)
{
  geocode <- validate_geocode(geocode)
  stopifnot(lubridate::is.POSIXct(time))
  stopifnot(length(tz) == 1)

  z <- list(nrow(geocode))
  for (i in 1:nrow(geocode)) {
    temp <- sun_angles_fast(time = time,
                            tz = tz,
                            geocode = dplyr::slice(geocode, i),
                            use.refraction = use.refraction)
    z[[i]] <- temp
  }
  # we supress warning of dropped attributes and restore them
  z <- suppressWarnings(dplyr::bind_rows(z))
  class(z[["solartime"]]) <- class(temp[["solartime"]])

  # first.iter <- TRUE
  # for (i in 1:nrow(geocode)) {
  #   temp <- sun_angles_fast(time = time,
  #                           tz = tz,
  #                           geocode = geocode[i, ],
  #                           use.refraction = use.refraction)
  #   if (first.iter) {
  #     z <- temp
  #     first.iter <- FALSE
  #   } else {
  #     z <- rbind(z, temp,
  #                make.row.names = FALSE)
  #   }
  # }

  # we use rbind instead of dplyr::bind_rows as the second drops the class attribute
  # for solartime.

  # assertion
  if (any(z[["elevation"]] < (-90)) || any(z[["elevation"]] > 90))
    warning("Returned 'elevation' value(s) off range")
  if (any(!is.na(z[["azimuth"]]) & (z[["azimuth"]] < -1e-10 |
         z[["azimuth"]] > (360 + 1e-10)))) {
      warning("Returned 'azimuth' values(s) off range")
  }
  z
}

#' @rdname sun_angles
#'
# Internal function, called by sun_angles()
#
sun_angles_fast <- function(time,
                            tz,
                            geocode,
                            use.refraction)
{
  # Input validation done in sun_angles() before calling this function.
  # stopifnot(!anyNA(time))
  # stopifnot(is.data.frame(geocode))
  # stopifnot(nrow(geocode == 1) && length(tz == 1))
  # We have a single geocode and all times are expressed in the same time zone!
  # If time is a vector we can vectorize the whole calculation, and do the
  # expensive calculations only once.

  lon <- geocode[["lon"]]

  lat <- geocode[["lat"]]

  address <- geocode[["address"]]

  cent <- julian_century(time)

  sun.lon.mean <- geom_mean_lon_sun(cent)
  sun.anom.mean <- geom_mean_anom_sun(cent)
  eccent.earth <- eccent_earth_orbit(cent)
  delta <- sun_eq_of_ctr(cent, sun.anom.mean)

  sun.lon <- sun.lon.mean + delta
#  sun.anom <- sun.anom.mean + delta
#  sun.dist <- sun_rad_vector(eccent.earth, sun.anom)
  sun.app.lon <- sun_app_lon(cent, sun.lon)
  sun.ecliptic <- mean_obliq_eclip(cent)
  obliq.corr <- obliq_corr(cent, sun.ecliptic)
#  rt.ascen <- sun_rt_ascen(sun.app.lon, obliq.corr)
  var.y <- var_y(obliq.corr)
  eq.of.time <- eq_of_time(mean.lon = sun.lon.mean,
                           eccent.earth = eccent.earth,
                           anom.mean = sun.anom.mean,
                           var.y = var.y)

  solar.time <- solar_tod(time, lat, lon, eq.of.time)
  hour.angle <- hour_angle(solar.time)
  sun.declin <- sun_decline(sun.app.lon, obliq.corr)
  zenith.angle <- zenith_angle(lat, hour.angle, sun.declin)
  elevation.angle <- 90 - zenith.angle
  if (use.refraction) {
    elevation.angle <-
    elevation.angle + atm_refraction_approx(elevation.angle)
  }
  azimuth.angle <- azimuth_angle(lat, hour.angle, zenith.angle, sun.declin)
  solar.time <- solar.time / 60 # hours
  class(solar.time) <- c("solar_time", class(solar.time))

  z <- tibble::tibble(time = lubridate::with_tz(time, tz),
                      tz = rep(tz, length(time)),
                      solartime = solar.time %% 24, # needed for DST
                      longitude = rep(lon, length(time)),
                      latitude = rep(lat, length(time)),
                      address = rep(address, length(time)),
                      azimuth = azimuth.angle,
                      elevation = elevation.angle,
                      declination = sun.declin,
                      eq.of.time = eq.of.time,
                      hour.angle = hour.angle,
                      .name_repair = "minimal")
  z
}

#' @rdname sun_angles
#'
#' @export
#'
sun_elevation <- function(time = lubridate::now(),
                          tz = lubridate::tz(time),
                          geocode = tibble::tibble(lon = 0,
                                                   lat = 51.5,
                                                   address = "Greenwich"),
                          use.refraction = FALSE)
{
  stopifnot(length(time) == 1 || nrow(geocode) == 1)
  sun_angles(time = time,
             tz = tz,
             geocode = geocode,
             use.refraction = use.refraction)[["elevation"]]
}

#' @rdname sun_angles
#'
#' @export
#'
sun_zenith_angle <- function(time = lubridate::now(),
                             tz = lubridate::tz(time),
                             geocode = tibble::tibble(lon = 0,
                                                      lat = 51.5,
                                                      address = "Greenwich"),
                             use.refraction = FALSE)
{
  stopifnot(length(time) == 1 || nrow(geocode) == 1)
  90 - sun_angles(time = time,
                  tz = tz,
                  geocode = geocode,
                  use.refraction = use.refraction)[["elevation"]]
}

#' @rdname sun_angles
#'
#' @export
#'
sun_azimuth <- function(time = lubridate::now(),
                        tz = lubridate::tz(time),
                        geocode = tibble::tibble(lon = 0,
                                                 lat = 51.5,
                                                 address = "Greenwich"),
                        use.refraction = FALSE)
{
  stopifnot(length(time) == 1 || nrow(geocode) == 1)
  sun_angles(time = time,
             tz = tz,
             geocode = geocode,
             use.refraction = use.refraction)[["azimuth"]]
}

#' Time difference between two time zones
#'
#' Returns the time difference in hours between two time zones at a given
#' instant in time.
#'
#' @param when datetime A time instant
#' @param tz.target,tz.reference character Two time zones using names
#' recognized by functions from package 'lubridate'
#'
#' @return A \code{numeric} value.
#'
#' @export
#'
tz_time_diff <- function(when = lubridate::now(),
                         tz.target = lubridate::tz(when),
                         tz.reference = "UTC") {
  if (lubridate::is.Date(when)) {
    when <- lubridate::as_datetime(when, tz = tz.target)
  }
  (as.numeric(lubridate::force_tz(when, tz.reference)) -
      as.numeric(lubridate::force_tz(when, tz.target))) / 3600
}

#' Times for sun positions
#'
#' Functions for calculating the timing of solar positions by means of function
#' \code{sun_angles}, given geographical coordinates and dates. They can be also
#' used to find the time for an arbitrary solar elevation between 90 and -90
#' degrees by supplying "twilight" angle(s) as argument.
#'
#' @param date "vector" of POSIXct times or Date objects, any valid TZ is allowed,
#'   default is current date at Greenwich.
#' @param tz character vector indicating time zone to be used in output.
#' @param geocode data frame with one or more rows and variables lon and lat as
#'   numeric values (degrees). If present, address will be copied to the output.
#' @param twilight character string, one of "none", "civil", "nautical",
#'   "astronomical", or a \code{numeric} vector of length one, or two, giving
#'   solar elevation angle(s) in degrees (negative if below the horizon).
#' @param unit.out character string, One of "datetime", "day", "hour", "minute",
#'   or "second".
#'
#' @return A tibble with variables day, tz, twilight.rise, twilight.set,
#'   longitude, latitude, address, sunrise, noon, sunset, daylength,
#'   nightlength or the corresponding individual vectors.
#'
#' @family astronomy related functions
#'
#' @details Twilight names are interpreted as follows. "none": solar elevation =
#'   0 degrees. "refraction": solar elevation = 0 degrees + refraction
#'   correction. "sunlight": upper rim of solar disk corrected for refraction,
#'   which the value used by the online NOAA Solar Calculator. "civil": -6
#'   degrees, "naval": -12 degrees, and "astronomical": -18 degrees. Unit names
#'   for output are as follows: "hours" times for sunrise and sunset are
#'   returned as times-of-day in hours since midnight. "date" or "datetime"
#'   return the same times as datetime objects with TZ set (this is much slower
#'   than the "hours"). Day length and night length are returned as numeric
#'   values expressed in hours when `"datetime"' is passed as argument to
#'   \code{unit.out}. If twilight is a numeric vector of length two, the element
#'   with index 1 is used for sunrise and that with index 2 for sunset.
#'
#' @seealso \code{\link{sun_angles}}.
#'
#' @note This function is an implementation of Meeus equations as used in NOAAs
#'   on-line web calculator, which are very precise and valid for a very broad
#'   range of dates. For sunrise and sunset the times are affected by refraction
#'   in the atmosphere, which does in turn depend on weather conditions. The
#'   effect of refraction on the apparent position of the sun is only an
#'   estimate based on "typical" conditions. The more tangential to the horizon
#'   is the path of the sun, the larger the effect of refraction is on the times
#'   of visual occlusion of the sun behind the horizon---i.e. the largest timing
#'   errors occur at high latitudes. The computation is not defined for
#'   latitudes 90 and -90 degrees, i.e. at the poles.
#'
#'   There exists a different R implementation of the same algorithms called
#'   "AstroCalcPureR" available as function \code{astrocalc4r} in package
#'   'fishmethods'. Although the equations used are almost all the same, the
#'   function signatures and which values are returned differ. In particular,
#'   the present implementation splits the calculation into two separate
#'   functions, one returning angles at given instants in time, and a separate
#'   one returning the timing of events for given dates.
#'
#'   In the current implementation functions \code{sunrise_time},
#'   \code{noon_time}, \code{sunset_time} and \code{day_length} are wrappers
#'   on \code{day_night}, so if more than one quantity is needed it is
#'   preferable to directly call \code{day_night} as it will be faster.
#'
#' @section Warning: Be aware that R's \code{Date} class does not save time zone
#'   metadata. This can lead to ambiguities in the current implementation as
#'   based on time instants. The argument passed to \code{date} should be
#'   of class \code{POSIXct}, in other words an instant in time, from which
#'   the correct date will be computed based on the \code{tz} argument.
#'
#' @export
#' @examples
#' library(lubridate)
#' my.geocode <- data.frame(lat = 60, lon = 25)
#' day_night(ymd("2015-05-30"), geocode = my.geocode)
#' day_night(ymd("2015-05-30") + days(1:10), geocode = my.geocode, twilight = "civil")
#' sunrise_time(ymd("2015-05-30"), geocode = my.geocode)
#' noon_time(ymd("2015-05-30"), geocode = my.geocode)
#' sunset_time(ymd("2015-05-30"), geocode = my.geocode)
#' day_length(ymd("2015-05-30"), geocode = my.geocode)
#' day_length(ymd("2015-05-30"), geocode = my.geocode, unit.out = "day")
#'
day_night <- function(date = lubridate::now(tzone = "UTC"),
                      tz = lubridate::tz(date),
                      geocode = tibble::tibble(lon = 0,
                                               lat = 51.5,
                                               address = "Greenwich"),
                      twilight = "none",
                      unit.out = "hours") {
  stopifnot(! anyNA(date))
  geocode <- validate_geocode(geocode)
  date <- as.Date(date, tz = "UTC")
#  date <- lubridate::floor_date(date, unit = "days") resulted in error!!

  if (unit.out == "date") {
    unit.out <- "datetime"
  } else if (unit.out == "days") {
    unit.out <- "day"
  } else if (unit.out == "hours") {
    unit.out <- "hour"
  } else if (unit.out == "minutes") {
    unit.out <- "minute"
  } else if (unit.out == "seconds") {
    unit.out <- "second"
  }

  z <- list(nrow(geocode))
  for (i in 1:nrow(geocode)) {
    temp <- day_night_fast(date = date,
                           tz = tz,
                           geocode = dplyr::slice(geocode, i),
                           twilight = twilight,
                           unit.out = unit.out)
    z[[i]] <- temp
  }
  # we supress warning of dropped attributes and restore them
  #  z <- suppressWarnings(dplyr::bind_rows(z))
  z <- do.call(rbind, z)

  # we use rbind instead of dplyr::bind_rows as the second drops the class attribute
  # for solartime.

  # assertion
  if (any(z[["daylength"]] < 0) || any(z[["nightlength"]] < 0)) {
    warning("Returned 'daylength/nightlength' value(s) off range")
  }
  attr(z, "unit.out") <- unit.out

  z
}

#' @rdname day_night
#'
day_night_fast <- function(date,
                           tz,
                           geocode,
                           twilight,
                           unit.out) {

  # Input validation done in day_night() before calling this function.
  # stopifnot(!anyNA(time))
  # stopifnot(is.data.frame(geocode))
  # stopifnot(nrow(geocode == 1) && length(tz == 1))
  # We have a single geocode and all dates are expressed in the same time zone!
  # As date is a "vector" we can vectorize the whole calculation, and do the
  # expensive calculations only once.

  tz <- tz[1]

  multiplier <- switch(unit.out,
                       hour = 1,
                       minute = 60,
                       second = 3600,
                       day = 1/24)

  # not vectorized, but possibly different angle for sunset and sunrise
  # twilight.angles is always of length 2
  twilight.angles <- twilight2angle(twilight)

  # geocode passed by day_night() has always one row
  if (!exists("address", geocode)) {
    geocode[["address"]] <- NA_character_
  }

  lon <- geocode[["lon"]]
  if (lon > 180 || lon < -180) {
    stop("Longitude is off-range.")
  }

  lat <- geocode[["lat"]]
  if (lat > 89.99 || lat < -89.99) {
    stop("Latitude is off-range.")
  }

  address <- geocode[["address"]]

#  date <- lubridate::floor_date(lubridate::with_tz(date, tzone = tz), unit = "day")
  date <- lubridate::as_date(date, tz = tz) # discards tz

   noon.of.date <- lubridate::as_datetime(date) + 43200 # as_datetime() is needed to obtain correct aswers!!
#   noon.of.date <- lubridate::as_datetime(date) + lubridate::seconds(43200) # faster
#  noon.of.date <- lubridate::with_tz(date, tzone = "UTC") + 43200 # faster

  cent <- julian_century(noon.of.date)

  tz.diff <- tz_time_diff(noon.of.date, tz.target = tz)

  sun.lon.mean <- geom_mean_lon_sun(cent)
  sun.anom.mean <- geom_mean_anom_sun(cent)
  eccent.earth <- eccent_earth_orbit(cent)
  delta <- sun_eq_of_ctr(cent, sun.anom.mean)

  sun.lon <- sun.lon.mean + delta
  sun.anom <- sun.anom.mean + delta
#  sun.dist <- sun_rad_vector(eccent.earth, sun.anom)
  sun.app.lon <- sun_app_lon(cent, sun.lon)
  sun.ecliptic <- mean_obliq_eclip(cent)
  obliq.corr <- obliq_corr(cent, sun.ecliptic)
#  rt.ascen <- sun_rt_ascen(sun.app.lon, obliq.corr)
  sun.declin <- sun_decline(sun.app.lon, obliq.corr)
  var.y <- var_y(obliq.corr)
  eq.of.time <- eq_of_time(mean.lon = sun.lon.mean,
                           eccent.earth = eccent.earth,
                           anom.mean = sun.anom.mean,
                           var.y = var.y)

  solar.noon <- solar_noon(lon, eq.of.time)

  # We need to test for 24 h and 0 h days
  sun.noon.elevation <- elevation_angle(lat, 0, sun.declin)
  low.sunrise.sun <- sun.noon.elevation < max(twilight.angles[1])
  low.sunset.sun <- sun.noon.elevation < max(twilight.angles[2])

  sun.midnight.elevation <- elevation_angle(lat, 180, sun.declin)
  high.sunrise.sun <- sun.midnight.elevation > min(twilight.angles[1])
  high.sunset.sun <- sun.midnight.elevation > min(twilight.angles[2])

  # Vectorized
  day.start <-
    ifelse(!low.sunrise.sun & !high.sunrise.sun,
           # normal case
           sunrise(solar.noon,
                   ha_sunrise(lat, sun.declin, nag = -twilight.angles[1])),
           0
    )
  sunrise <- ifelse(day.start == 0, NA_real_, day.start)

  day.end <-
    ifelse(!low.sunset.sun & !high.sunset.sun,
           # normal case
           sunset(solar.noon,
                  ha_sunrise(lat, sun.declin, nag = -twilight.angles[2])),
           1
    )
  sunset <- ifelse(day.end == 1, NA_real_, day.end)

  daylength.hours <- ifelse(low.sunrise.sun & low.sunset.sun,
                            0,
                            (day.end - day.start) * 24)

  # we assemble the data frame for one geographic location
  # the slowest part of calculations are the calls to lubridate
  # so we try to minimize them according to output format.

  if (unit.out == "datetime") {
    sunrise.time <- lubridate::as_datetime(date, tz = tz) +
      lubridate::seconds(sunrise * 86400)
    noon.time    <- lubridate::as_datetime(date, tz = tz) +
      lubridate::seconds(solar.noon * 86400)
    sunset.time  <- lubridate::as_datetime(date, tz = tz) +
      lubridate::seconds(sunset * 86400)

    z <- tibble::tibble(day           = date,
                        tz            = rep(tz, length(date)),
                        twilight.rise = rep(twilight.angles[1], length(date)),
                        twilight.set  = rep(twilight.angles[2], length(date)),
                        longitude     = rep(lon, length(date)),
                        latitude      = rep(lat, length(date)),
                        address       = rep(address, length(date)),
                        sunrise       = sunrise.time, #lubridate::with_tz(sunrise.time, tzone = tz),
                        noon          = noon.time, #lubridate::with_tz(noon.time, tzone = tz),
                        sunset        = sunset.time, #lubridate::with_tz(sunset.time, tzone = tz),
                        daylength     = daylength.hours,
                        nightlength   = 24 - daylength.hours,
                        .name_repair  = "minimal"
    )
  } else if (unit.out %in% c("day", "hour", "minute", "second")) {
    sunrise.tod <- (sunrise * 24 + tz.diff) %% 24
    noon.tod <- (solar.noon * 24 + tz.diff) %% 24
    sunset.tod <- (sunset * 24 + tz.diff) %% 24

    z <- tibble::tibble(day           = date,
                        tz            = rep(tz, length(date)),
                        twilight.rise = rep(twilight.angles[1], length(date)),
                        twilight.set  = rep(twilight.angles[2], length(date)),
                        longitude     = rep(lon, length(date)),
                        latitude      = rep(lat, length(date)),
                        address       = rep(address, length(date)),
                        sunrise       = sunrise.tod * multiplier,
                        noon          = noon.tod * multiplier,
                        sunset        = sunset.tod * multiplier,
                        daylength     = daylength.hours * multiplier,
                        nightlength   = (24 - daylength.hours) * multiplier,
                        .name_repair  = "minimal"
    )
  } else {
    stop("Unit out '", unit.out, "' not recognized")
  }

  z
}

#' twilight argument check and conversion
#'
#' @return numeric  Solar elevation angle at sunrise or sunset
#' @keywords internal
twilight2angle <- function(twilight) {
  if (!is.numeric(twilight)) {
    if (twilight == "none") { # center of solar disk
      twilight_angle <- c(0, 0)
    } else if (twilight == "refraction") { #  center of solar disk, refraction corrected
      twilight_angle <- c(-0.4819444, -0.4819444)
    } else if (twilight == "sunlight") { # upper rim of solar disk, refraction corrected
      twilight_angle <- c(-0.833, -0.833)
    } else if (twilight == "civil") { # refraction corrected
      twilight_angle <- c(-6, -6)
    } else if (twilight == "nautical") { # refraction corrected
      twilight_angle <- c(-12, -12)
    } else if (twilight == "astronomical") { # refraction corrected
      twilight_angle <- c(-18, -18)
    } else {
      twilight_angle <- c(NA, NA)
    }
  } else {
    if (length(twilight) == 1) {
      twilight_angle <- rep(twilight, 2)
    } else if (length(twilight) == 2) {
      twilight_angle <- twilight
    } else {
      twilight_angle <- c(NA, NA)
    }
    twilight_angle <- ifelse(twilight_angle < 90, twilight_angle, NA)
    twilight_angle <- ifelse(twilight_angle > -90, twilight_angle, NA)
  }
  if (any(is.na(twilight_angle))) {
    stop("Unrecognized argument value for 'twilight': ", twilight)
  }
  twilight_angle
}

#' @rdname day_night
#' @export
#' @return \code{noon_time}, \code{sunrise_time} and \code{sunset_time} return a
#'   vector of POSIXct times
noon_time <- function(date = lubridate::today(),
                      tz = lubridate::tz(date),
                      geocode = tibble::tibble(lon = 0,
                                               lat = 51.5,
                                               address = "Greenwich"),
                      twilight = "none",
                      unit.out = "datetime") {
  stopifnot(length(date) == 1 || nrow(geocode) == 1)
  day_night(date = date,
            tz = tz,
            geocode = geocode,
            twilight = twilight,
            unit.out = unit.out)[["noon"]]
}

#' @rdname day_night
#'
#' @export
sunrise_time <- function(date = lubridate::today(),
                         tz = lubridate::tz(date),
                         geocode = tibble::tibble(lon = 0,
                                                  lat = 51.5,
                                                  address = "Greenwich"),
                         twilight = "sunlight", unit.out = "datetime") {
 #  stopifnot(length(date) == 1L || nrow(geocode) == 1L)
  day_night(date = date,
            tz = tz,
            geocode = geocode,
            twilight = twilight,
            unit.out = unit.out)[["sunrise"]]
}

#' @rdname day_night
#' @export
#'
sunset_time <- function(date = lubridate::today(),
                        tz = lubridate::tz(date),
                        geocode = tibble::tibble(lon = 0,
                                                 lat = 51.5,
                                                 address = "Greenwich"),
                        twilight = "sunlight", unit.out = "datetime") {
  # stopifnot(length(date) == 1L || nrow(geocode) == 1L)
  day_night(date = date,
            tz = tz,
            geocode = geocode,
            twilight = twilight,
            unit.out = unit.out)[["sunset"]]
}

#' @rdname day_night
#'
#' @export
#' @return \code{day_length} and \code{night_length} return numeric a vector
#'   giving the length in hours
day_length <- function(date = lubridate::now(),
                       tz = "UTC",
                       geocode = tibble::tibble(lon = 0,
                                                lat = 51.5,
                                                address = "Greenwich"),
                       twilight = "sunlight", unit.out = "hours") {
  # stopifnot(length(date) == 1L || nrow(geocode) == 1L)
  day_night(date = date,
            tz = tz,
            geocode = geocode,
            twilight = twilight,
            unit.out = unit.out)[["daylength"]]
}

#' @rdname day_night
#'
#' @export
#' @note \code{night_length} returns the length of night-time conditions in one
#'   day (00:00:00 to 23:59:59), rather than the length of the night between two
#'   consecutive days.
night_length <- function(date = lubridate::now(),
                         tz = "UTC",
                         geocode = tibble::tibble(lon = 0,
                                                  lat = 51.5,
                                                  address = "Greenwich"),
                         twilight = "sunlight", unit.out = "hours") {
  # stopifnot(length(date) == 1L || nrow(geocode) == 1L)
  day_night(date = date,
              tz = tz,
              geocode = geocode,
              twilight = twilight,
              unit.out = unit.out)[["nightlength"]]
}

#' Convert date to time-of-day in hours, minutes or seconds
#'
#' @param x a datetime object accepted by lubridate functions
#' @param unit.out character string, One of "datetime", "hour", "minute", or "second".
#' @param tz character string indicating time zone to be used in output.
#'
#' @export
#'
as_tod <- function(x, unit.out = "hours", tz = NULL) {
  stopifnot(lubridate::is.timepoint(x))
  if (!is.null(tz)) {
    x <- lubridate::with_tz(x, tzone = tz[1])
  }
  if (unit.out == "hours") {
    lubridate::hour(x) + lubridate::minute(x) / 60 + lubridate::second(x) / 3600
  } else if (unit.out == "minutes") {
    lubridate::hour(x) * 60 + lubridate::minute(x) + lubridate::second(x) / 60
  } else if (unit.out == "seconds") {
    lubridate::hour(x) * 3600 + lubridate::minute(x) * 60 + lubridate::second(x)
  } else {
    stop("Unrecognized 'unit.out': ", unit.out)
  }
}

#' Local solar time
#'
#' \code{solar_time} computes from a time and geocode, the time of day expressed
#' in seconds since midnight. \code{solar_date} returns the same instant in time
#' as a date-time object. Solar time is useful when we want to plot data
#' according to the local solar time of day, irrespective of the date. Solar
#' date is useful when we want to plot a time series stretching for several days
#' using the local solar time but distinguishing between days.
#'
#' @param time POSIXct Time, any valid time zone (TZ) is allowed, default is
#'   current time
#' @param geocode data frame with variables lon and lat as numeric values
#'   (degrees).
#' @param unit.out character string, One of "datetime", "hour", "minute", or "second".
#'
#' @family astronomy related functions
#'
#' @note The algorithm is approximate, it calculates the difference between
#'   local solar noon and noon in the time zone of \code{time} and uses this
#'   value for the whole day when converting times into solar time. Days are not
#'   exactly 24 h long. Between successive days the shift is only a few seconds,
#'   and this leads to a small jump at midnight.
#'
#' @section Warning!:
#'   Returned values are computed based on the time zone of the argument for
#'   parameter time. In the case of solar time, this timezone does not affect
#'   the result. However, in the case of solar dates the date part may be
#'   off by one day, if the time zone does not match the coordinates of the
#'   geocode value provided as argument.
#'
#' @return For \code{solar_time()} numeric value in seconds from midnight but
#'   with an additional class attribute "solar.time".
#'
#' @export
#'
#' @examples
#' BA.geocode <-
#'   data.frame(lon = -58.38156, lat = -34.60368, address = "Buenos Aires, Argentina")
#' sol_t <- solar_time(lubridate::dmy_hms("21/06/2016 10:00:00", tz = "UTC"),
#'                     BA.geocode)
#' sol_t
#' class(sol_t)
#'
#' sol_d <- solar_time(lubridate::dmy_hms("21/06/2016 10:00:00", tz = "UTC"),
#'                     BA.geocode,
#'                     unit.out = "datetime")
#' sol_d
#' class(sol_d)
#'
solar_time <- function(time = lubridate::now(),
                       geocode = tibble::tibble(lon = 0,
                                                lat = 51.5,
                                                address = "Greenwich"),
                       unit.out = "time")
{
  # solar time in hours from midnight
  solar.time <- sun_angles(time = time,
                           tz = lubridate::tz(time),
                           geocode = geocode)[["solartime"]]
  switch(unit.out,
         "date" = as.solar_date(solar.time, time),
         "datetime" = as.solar_date(solar.time, time),
         "time" = solar.time,
         "days" = as.numeric(solar.time) / 24,
         "hours" = as.numeric(solar.time),
         "minutes" = as.numeric(solar.time * 60),
         "seconds" = as.numeric(solar.time * 3600)
  )
}

#' Convert a solar_time object into solar_date object
#'
#' @param x solar_time object.
#' @param time an R date time object
#'
#' @return For method \code{as.solar_date()} a date-time object with the class attr
#'   set to "solar.time". This is needed only for unambiguous formatting and
#'   printing.
#'
#' @export
#'
as.solar_date <- function(x, time)
{
  stopifnot(is.solar_time(x))
  stopifnot(lubridate::is.timepoint(time))
  solar.date <-
    lubridate::floor_date(time, unit = "days") +
    lubridate::seconds(as.numeric(x) * 3600)
  class(solar.date) <- c("solar_date", class(solar.date))
  solar.date
}

#' Query class
#'
#' @param x an R object.
#'
#' @family astronomy related functions
#'
#' @export
#'
is.solar_time <- function(x) {
  inherits(x, "solar_time")
}

#' @rdname is.solar_time
#'
#' @export
#'
is.solar_date <- function(x) {
  inherits(x, "solar_date")
}

#' Encode in a Common Format
#'
#' Format a \code{solar_time} object for pretty printing
#'
#' @param x an R object
#' @param ... ignored
#' @param sep character used as separator
#'
#' @family astronomy related functions
#'
#' @export
#'
format.solar_time <- function(x, ..., sep = ":") {
  hours <- as.integer(trunc(x))
  minutes <- as.integer((x * 60) %% 60)
  seconds <- as.integer((x * 3600) %% 60)
  fmt <- paste(rep("%02d", 3), collapse = sep)
  time_string <-
    sprintf(fmt = fmt, hours, minutes, seconds)
  time_string
}

#' Print solar time and solar date objects
#'
#' @param x an R object
#' @param ... passed to \code{format} method
#'
#' @family astronomy related functions
#'
#' @note Default is to print the underlying POSIXct as a solar time.
#'
#' @export
#'
print.solar_time <- function(x, ...) {
  print(format(x, ...))
  invisible(x)
}

#' @rdname print.solar_time
#'
#' @export
#'
print.solar_date <- function(x, ...) {
  print(paste(format(x, ...), "solar"))
  invisible(x)
}

#' Validate a geocode
#'
#' Convert to tibble, check data bounds, convert address to character if
#' it is not, or add character NAs if the address column is missing.
#'
#' @keywords internal
#'
validate_geocode <- function(geocode) {
  geocode <- tibble::as_tibble(geocode, .name_repair = "minimal")
  stopifnot(nrow(geocode) >= 1) # needs to be replace by generation of no output in all fucntions
  if (any(geocode[["lon"]] > 180 | geocode[["lon"]] < -180)) {
    stop("Longitude is off-range.")
  }
  if (any(geocode[["lat"]] > 89.99 | geocode[["lat"]] < -89.99)) {
    stop("Latitude is off-range.")
  }
  if (!exists("address", geocode)) {
    geocode[["address"]] <- NA_character_
  } else if (!is.character(geocode[["address"]])) {
    geocode[["address"]] <- as.character(geocode[["address"]])
  }
  geocode
}
