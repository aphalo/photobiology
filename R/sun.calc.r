#' Solar angles
#'
#' This function returns the solar angles for a given time and location.
#'
#' @param time POSIXct Time, any valid time zone (TZ) is allowed, default is
#'   current time
#' @param tz character string indicating time zone to be used in output.
#' @param geocode data frame with variables lon and lat as numeric values
#'   (degrees).
#' @param use.refraction logical Flag indicating whether to correct for
#'   fraction in the atmosphere
#'
#' @return A data.frame with variables time (in same TZ as input), TZ, solartime,
#'   longitude, latitude, address, azimuth, and elevation.
#'
#' @family astronomy related functions
#'
#' @export
#' @examples
#' library(lubridate)
#' sun_angles()
#' sun_angles(ymd_hms("2014-09-23 12:00:00"))
#' sun_angles(ymd_hms("2014-09-23 12:00:00"),
#'            geocode = data.frame(lat=60, lon=0))
#'
sun_angles <- function(time = lubridate::now(),
                       tz = lubridate::tz(time),
                       geocode = data.frame(lon = 0, lat = 51.5, address = "Greenwich"),
                       use.refraction = FALSE)
{
  stopifnot(!anyNA(time))
  stopifnot(is.data.frame(geocode))

  # if time is a vector or list for convenience we vectorize
  if (length(time) > 1) {
    first.iter <- TRUE
    for (i in 1:length(time)) {
      zz <- sun_angles(time = time[i],
                      tz = tz,
                      geocode = geocode,
                      use.refraction = use.refraction)
      if (first.iter) {
        z <- zz
        first.iter <- FALSE
      } else {
        z <- rbind(z, zz)
      }
    }
    return(z)
  }

  # from here onwards we are dealing with a single date but possibly
  # several time zones, so for the time being we repeat all calculations
  # for each location

  if (length(tz) == 1 && nrow(geocode) > 1) {
    tz <- rep(tz, nrow(geocode))
  } else if (length(tz) != nrow(geocode)) {
    stop("'tz' argument of wrong length")
  }

  # vectorized over geocodes

  if (!exists("address", geocode)) {
    geocode[["address"]] <- NA_character_
  }

  first.iter <- TRUE
  for (i in 1:nrow(geocode)) {
    lon <- geocode[i, "lon"]
    lat <- geocode[i, "lat"]
    address <- geocode[i, "address"]

    if (i == 1 || previous.tz != tz[i]) {

      previous.tz <- tz[i]

      # only needs recalculation when tz changes

      cent <- julian_century(time)

      sun.lon.mean <- geom_mean_lon_sun(cent)
      sun.anom.mean <- geom_mean_anom_sun(cent)
      eccent.earth <- eccent_earth_orbit(cent)
      delta <- sun_eq_of_ctr(cent, sun.anom.mean)

      sun.lon <- sun.lon.mean + delta
      sun.anom <- sun.anom.mean + delta
      sun.dist <- sun_rad_vector(eccent.earth, sun.anom)
      sun.app.lon <- sun_app_lon(cent, sun.lon)
      sun.ecliptic <- mean_obliq_eclip(cent)
      obliq.corr <- obliq_corr(cent, sun.ecliptic)
      rt.ascen <- sun_rt_ascen(sun.app.lon, obliq.corr)
      sun.declin <- sun_declin(sun.app.lon, obliq.corr)
      var.y <- var_y(obliq.corr)
      eq.of.time <- eq_of_time(mean.lon = sun.lon.mean,
                               eccent.earth = eccent.earth,
                               anom.mean = sun.anom.mean,
                               var.y = var.y)

    }

    solar.time <- solar_tod(time, lat, lon, eq.of.time)
    hour.angle <- hour_angle(solar.time)
    zenith.angle <- zenith_angle(lat, hour.angle, sun.declin)
    elevation.angle <- 90 - zenith.angle
    if (use.refraction) {
      elevation.angle <-
        elevation.angle + atm_refraction_approx(elevation.angle)
    }
    azimuth.angle <- azimuth_angle(lat, hour.angle, zenith.angle, sun.declin)

    solar.time <- solar.time / 60 # hours
    class(solar.time) <- c("solar_time", class(solar.time))

    yy <- data.frame(time = lubridate::with_tz(time, tz[i]),
                     tz = tz[i],
                     solartime = solar.time,
                     longitude = lon,
                     latitude = lat,
                     address = address,
                     azimuth = azimuth.angle,
                     elevation = elevation.angle)

    # we bind the data frames together
    if (first.iter) {
      y <- yy
      first.iter <- FALSE
    } else {
      y <- rbind(y, yy)
    }
  }
  # assertion

  if (any(y[["elevation"]] < (-90)) || any(y[["elevation"]] > 90))
    stop("Returned 'elevation' value(s) off range")
  if (any(y[["azimuth"]] < 0) || any(y[["azimuth"]] > 360))
    stop("Returned 'azimuth' values(s) off range")

  y
}

#' @rdname sun_angles
#'
#' @export
#'
sun_elevation <- function(time = lubridate::now(),
                          tz = lubridate::tz(time),
                          geocode = data.frame(lon = 0, lat = 51.5, address = "Greenwich"),
                          use.refraction = FALSE)
{
  sun_angles(time = time,
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
                          geocode = data.frame(lon = 0, lat = 51.5, address = "Greenwich"),
                          use.refraction = FALSE)
{
  sun_angles(time = time,
             tz = tz,
             geocode = geocode,
             use.refraction = use.refraction)[["elevation"]] - 90
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
#' @export
#'
tz_time_diff <- function(when = lubridate::now(),
                         tz.target = Sys.timezone(),
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
#' @param date vector of POSIXct times or Date objects, any valid TZ is allowed,
#'   default is current date
#' @param tz vector of character string indicating time zone to be used in output.
#' @param geocode data frame with one or more rows and variables lon and lat as
#'   numeric values (degrees).
#' @param twilight character string, one of "none", "civil", "nautical",
#'   "astronomical", or a \code{numeric} vector of length one, or two, giving
#'   solar elevation angle(s) in degrees (negative if below the horizon).
#' @param unit.out charater string, One of "datetime", "hour", "minute", or "second".
#'
#' @return A data.fraame with variables day, tz, twilight.rise, twilight.set,
#'   longitude, latitude, address, sunrise, noon, sunset, daylength,
#'   nightlength.
#'
#' @note If twilight is a numeric vector of length two, the element with index 1
#'   is used for sunrise and that with index 2 for sunset.
#'
#' @family astronomy related functions
#'
#' @details Twilight names are interpreted as follows. "none": solar elevation =
#'   0 degrees. "refraction": solar elevation = 0 degrees + refraction
#'   correction. "sunlight": upper rim of solar disk corrected for refraction.
#'   "civil": -6 degrees, "naval": -12 degrees, and "astronomical": -18 degrees.
#'   Unit names for output are as follows: "hours" times for sunrise and sunset
#'   are returned as times-of-day in hours since middnight. "date" or "datetime"
#'   return the same times as datetime objects with TZ set (this is much slower
#'   the "hours"). Day length and night length are always returned as numeric
#'   values expressed in hours.
#'
#' @export
#' @examples
#' library(lubridate)
#' my.geocode <- data.frame(lat = 60, lon = 25)
#' day_night(ymd("2015-05-30"), geocode = my.geocode, twilight = "civil")
#'
day_night <- function(date = lubridate::today(),
                      tz = Sys.timezone(),
                      geocode = data.frame(lon = 0, lat = 51.5, address = "Greenwich"),
                      twilight = "none",
                      unit.out = "hours") {
  stopifnot(! anyNA(date))
  stopifnot(is.data.frame(geocode))
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

  multiplier <- switch(unit.out,
                       hour = 1,
                       minute = 60,
                       second = 3600,
                       day = 1/24)

  # if date is a vector or list for convenience we vectorize
  if (length(date) > 1) {
    first.iter <- TRUE
    for (i in 1:length(date)) { # using (d in date) d is not a Date object!
      zz <- day_night(date = date[i],
                      tz = tz,
                      geocode = geocode,
                      twilight = twilight,
                      unit.out = unit.out)
      if (first.iter) {
        z <- zz
        first.iter <- FALSE
      } else {
        z <- rbind(z, zz)
      }
    }
    return(z)
  }

  # from here onwards we are dealing with a single date but possibly
  # several time zones, so for the time being we repeat all calculations
  # for each location

  if (unit.out == "datetime") {
    duration.unit.out <- "hours"
  } else {
    duration.unit.out <- unit.out
  }

  if (length(tz) == 1 && nrow(geocode) > 1) {
    tz <- rep(tz, nrow(geocode))
  } else if (length(tz) != nrow(geocode)) {
    stop("'tz' argument of wrong length")
  }

  # not vectorized
  twilight.angles <- twilight2angle(twilight)

  # vectorized over geocodes

  if (!exists("address", geocode)) {
    geocode[["address"]] <- NA_character_
  }

  first.iter <- TRUE
  for (i in 1:nrow(geocode)) {
    lon <- geocode[i, "lon"]
    lat <- geocode[i, "lat"]
    address <- geocode[i, "address"]

    if (i == 1 || previous.tz != tz[i]) {

      previous.tz <- tz[i]

      # only needs recalculation when tz changes

      date <- lubridate::as_date(date, tz = tz[i])

      noon.of.date <- lubridate::as_datetime(date) + lubridate::hours(12)
      cent <- julian_century(noon.of.date)

      tz.diff <- tz_time_diff(noon.of.date, tz.target = tz[i])

      sun.lon.mean <- geom_mean_lon_sun(cent)
      sun.anom.mean <- geom_mean_anom_sun(cent)
      eccent.earth <- eccent_earth_orbit(cent)
      delta <- sun_eq_of_ctr(cent, sun.anom.mean)

      sun.lon <- sun.lon.mean + delta
      sun.anom <- sun.anom.mean + delta
      sun.dist <- sun_rad_vector(eccent.earth, sun.anom)
      sun.app.lon <- sun_app_lon(cent, sun.lon)
      sun.ecliptic <- mean_obliq_eclip(cent)
      obliq.corr <- obliq_corr(cent, sun.ecliptic)
      rt.ascen <- sun_rt_ascen(sun.app.lon, obliq.corr)
      sun.declin <- sun_declin(sun.app.lon, obliq.corr)
      var.y <- var_y(obliq.corr)
      eq.of.time <- eq_of_time(mean.lon = sun.lon.mean,
                               eccent.earth = eccent.earth,
                               anom.mean = sun.anom.mean,
                               var.y = var.y)

    }

    solar.noon <- solar_noon(lon, eq.of.time)

    # We need to test for 24 h and 0 h days
    sun.noon.elevation <- elevation_angle(lat, 0, sun.declin)
    sun.midnight.elevation <- elevation_angle(lat, 180, sun.declin)

    if (sun.noon.elevation < max(twilight.angles)) {
      # at least one
      low.sun <- TRUE
      high.sun <- FALSE
      if (sun.noon.elevation < twilight.angles[1]) {
        sunrise <- NA_real_
      } else {
        ha.sunrise1 <- ha_sunrise(lat, sun.declin, ang = -twilight.angles[1])
        sunrise <- sunrise(solar.noon, ha.sunrise1)
      }
      if (sun.noon.elevation < twilight.angles[2]) {
        sunset <- NA_real_
      } else {
        ha.sunrise2 <- ha_sunrise(lat, sun.declin, ang = -twilight.angles[2])
        sunset <- sunset(solar.noon, ha.sunrise2)
      }
    } else if (sun.midnight.elevation > min(twilight.angles)) {
      # at least one
      low.sun <- FALSE
      high.sun <- TRUE
      if (sun.midnight.elevation > twilight.angles[1]) {
        sunrise <- NA_real_
      } else {
        ha.sunrise1 <- ha_sunrise(lat, sun.declin, ang = -twilight.angles[1])
        sunrise <- sunrise(solar.noon, ha.sunrise1)
      }
      if (sun.midnight.elevation > twilight.angles[2]) {
        sunset <- NA_real_
      } else {
        ha.sunrise2 <- ha_sunrise(lat, sun.declin, ang = -twilight.angles[2])
        sunset <- sunset(solar.noon, ha.sunrise2)
      }
    } else if (twilight.angles[1] == twilight.angles[2]) {
      low.sun <- FALSE
      high.sun <- FALSE
      ha.sunrise <- ha_sunrise(lat, sun.declin, ang = -twilight.angles[1])
      sunrise <- sunrise(solar.noon, ha.sunrise)
      sunset <- sunset(solar.noon, ha.sunrise)
    } else {
      low.sun <- FALSE
      high.sun <- FALSE
      ha.sunrise1 <- ha_sunrise(lat, sun.declin, ang = -twilight.angles[1])
      ha.sunrise2 <- ha_sunrise(lat, sun.declin, ang = -twilight.angles[2])
      sunrise <- sunrise(solar.noon, ha.sunrise1)
      sunset <- sunset(solar.noon, ha.sunrise2)
    }

    # we assemble the data frame for one geographic location
    # the slowest part of calcualtions are the calls to lubridate
    # so we try to minize them according to output format.

    if (unit.out == "datetime") {
      sunrise.time <- lubridate::as_datetime(date, tz = tz[i]) +
        lubridate::seconds(sunrise * 86400)
      noon.time    <- lubridate::as_datetime(date, tz = tz[i]) +
        lubridate::seconds(solar.noon * 86400)
      sunset.time  <- lubridate::as_datetime(date, tz = tz[i]) +
        lubridate::seconds(sunset * 86400)

      if (low.sun && (is.na(sunset.time) || is.na(sunset.time))) {
        daylength.hours <- 0
      } else if (high.sun && is.na(sunset.time) && is.na(sunset.time)) {
        daylength.hours <- 24
      } else if (high.sun) {
        if (is.na(sunrise.time)) {
          daylength.duration <- sunset.time - lubridate::as_datetime(date, tz = tz[i])
        } else if (is.na(sunset.time)) {
          daylength.duration <- lubridate::as_datetime(date, tz = tz[i]) +
            lubridate::seconds(24 * 86400) - sunrise.time
        }
        daylength.hours <- as.numeric(daylength.duration, duration.unit.out)
      } else {
        daylength.duration <- sunset.time - sunrise.time
        daylength.hours <- as.numeric(daylength.duration, duration.unit.out)
      }


      yy <- data.frame(day           = date,
                       tz            = tz[i],
                       twilight.rise = twilight[1],
                       twilight.set  = ifelse(length(twilight) == 1,
                                              twilight[1], twilight[2]),
                       longitude     = lon,
                       latitude      = lat,
                       address       = address,
                       sunrise       = lubridate::with_tz(sunrise.time, tzone = tz[i]),
                       noon          = lubridate::with_tz(noon.time, tzone = tz[i]),
                       sunset        = lubridate::with_tz(sunset.time, tzone = tz[i]),
                       daylength     = daylength.hours,
                       nightlength   = 24 - daylength.hours
      )
    } else if (unit.out %in% c("day", "hour", "minute", "second")) {
      sunrise.tod <- (sunrise * 24 + tz.diff) %% 24
      noon.tod <- (solar.noon * 24 + tz.diff) %% 24
      sunset.tod <- (sunset * 24 + tz.diff) %% 24

      if (low.sun && (is.na(sunset.tod) || is.na(sunset.tod))) {
        daylength.hours <- 0
      } else if (high.sun && is.na(sunset.tod) && is.na(sunset.tod)) {
        daylength.hours <- 24
      } else if (high.sun) {
        if (is.na(sunrise.tod)) {
          daylength.hours <- sunset.tod
        } else if (is.na(sunset.tod)) {
          daylength.hours <- 24 - sunrise.tod
        }
      } else {
        daylength.hours <- sunset.tod - sunrise.tod
      }

      yy <- data.frame(day           = date,
                       tz            = tz[i],
                       twilight.rise = twilight[1],
                       twilight.set  = ifelse(length(twilight) == 1,
                                              twilight[1], twilight[2]),
                       longitude     = lon,
                       latitude      = lat,
                       address       = address,
                       sunrise       = sunrise.tod * multiplier,
                       noon          = noon.tod * multiplier,
                       sunset        = sunset.tod * multiplier,
                       daylength     = daylength.hours * multiplier,
                       nightlength   = (24 - daylength.hours) * multiplier
      )
    } else {
      stop("Unit out '", unit.out, "' not recognized")
    }

    # we bind the data frames together
    if (first.iter) {
      y <- yy
      first.iter <- FALSE
    } else {
      y <- rbind(y, yy)
    }
  }
  attr(y, "unit.out") <- unit.out
  y
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
                      tz = Sys.timezone(),
                      geocode = data.frame(lon = 0, lat = 51.5, address = "Greenwich"),
                      twilight = "none",
                      unit.out = "datetime") {
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
                         tz = Sys.timezone(),
                         geocode = data.frame(lon = 0, lat = 51.5, address = "Greenwich"),
                         twilight = "sunlight", unit.out = "datetime") {
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
                        tz = Sys.timezone(),
                        geocode = data.frame(lon = 0, lat = 51.5, address = "Greenwich"),
                        twilight = "sunlight", unit.out = "datetime") {
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
day_length <- function(date = lubridate::today(),
                       tz = "UTC",
                       geocode = data.frame(lon = 0, lat = 51.5, address = "Greenwich"),
                       twilight = "sunlight", unit.out = "hours") {
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
#'   consequtive days.
night_length <- function(date = lubridate::today(),
                         tz = "UTC",
                         geocode = data.frame(lon = 0, lat = 51.5, address = "Greenwich"),
                         twilight = "sunlight", unit.out = "hours") {
  day_night(date = date,
              tz = tz,
              geocode = geocode,
              twilight = twilight,
              unit.out = unit.out)[["nightlength"]]
}

#' Convert date to time-of-day in hours, minutes or seconds
#'
#' @param x a datetime object accepted by lubridate functions
#' @param unit.out charater string, One of "datetime", "hour", "minute", or "second".
#' @param tz character string indicating time zone to be used in output.
#'
#' @export
#'
as_tod <- function(x, unit.out = "hours", tz = NULL) {
  stopifnot(lubridate::is.timepoint(x))
  if (!is.null(tz)) {
    x <- lubridate::with_tz(x, tzone = tz)
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
#' date is useful when we want to plot a time series streching for several days
#' using the local solar time but distinguishing between days.
#'
#' @param time POSIXct Time, any valid time zone (TZ) is allowed, default is
#'   current time
#' @param geocode data frame with variables lon and lat as numeric values
#'   (degrees).
#' @param unit.out charater string, One of "datetime", "hour", "minute", or "second".
#'
#' @family astronomy related functions
#'
#' @note The algorithm is approximate, it calculates the difference between
#'   local solar noon and noon in the time zone of \code{time} and uses this
#'   value for the whole day when converting times into solar time. Days are not
#'   exactly 24 h long. Between succesive days the shift is only a few seconds,
#'   and this leads to a small jump at midnight.
#'
#' @section Warning!:
#'   Returned values are computed based on the time zone of the argument for
#'   parameter time. In the case of solar time, this timezone does not affect
#'   the result. However, in the case of solar dates the date part may be be
#'   off by one day, if the time zone does not match the coordinates of the
#'   goecode value provided as argument.
#'
#' @return For \code{solar_time()} numeric value in seconds from midnight but
#'   with an additional class attribute "solar.time".
#'
#' @export
#'
#' @examples
#' # BA.geocode <- ggmap::geocode("Buenos Aires, Argentina")
#' BA.geocode <- data.frame(lon = -58.38156, lat = -34.60368)
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
                       geocode = data.frame(lon = 0, lat = 51.5, address = "Greenwich"),
                       unit.out = "time")
{
  # solar time in hours from midnight
  solar.time <- sun_angles(time = time,
                           tz = Sys.timezone(),
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
#'   set to "solar.time". This is needed only for unambiguous formating and
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
