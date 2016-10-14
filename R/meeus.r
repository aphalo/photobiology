#' Solar astronomy using Meeus' algorithm
#'
#' Low level functions based on NOAA's Excell worksheet
#'
#' @param time dateTime
#' @param x numeric Julian century
#' @param anom numeric Solar anomaly in degrees
#' @param eccent numric Eccentricity of Earth orbit
#' @param eclip numeric Ecliptic
#' @param app.lon,obliq.corr,mean.lon,ang numeric
#' @param eq.of.time,ha.sunrise,noon numeric
#' @param lat,lon numeric Geographic coordinates in degrees
#'
#' @keywords internal
#'
#' @examples
#' format(julian.day(dmy_hms("01/01/2010 12:00:00")), digits = 21)
#' cent <- julian.century(dmy_hms("01/01/2010 12:00:00"))
#'
#' sun.lon.mean <- geom_mean_lon_sun(cent)
#' sun.anom.mean <- geom_mean_anom_sun(cent)
#' eccent.earth <- eccent_earth_orbit(cent)
#' delta <- sun_eq_of_ctr(cent, anom.mean)
#'
#' sun.lon <- sun.lon.mean + delta
#' sun.anom <- sun.anom.mean + delta
#' sun.dist <- sun_rad_vector(eccent.earth, sun.anom)
#' sun.app.lon <- sun_app_lon(cent, sun.lon)
#' sun.ecliptic <- mean_obliq_eclip(cent)
#' obliq.corr <- obliq_corr(cent, sun.ecliptic)
#' rt.ascen <- sun_rt_ascen(sun.app.lon, obliq.corr)
#' sun.declin <- sun_declin(sun.app.lon, obliq.corr)
#' var.y <- var_y(obliq.corr)
#' eq.of.time <- eq_of_time(sun.lon.mean,
#'                          eccent.earth,
#'                          sun.anom.mean,
#'                          var.y)
#' ha.sunrise <- ha_sunrise(40, sun.declin)
#' solar.noon.utc <- solar_noon_utc(-105, eq.of.time)
#' sunrise.utc <- sunrise_utc(solar.noon.utc)
#' sunlight.duration <- sunlight_duration(ha.sunrise)
#' true_solar_time(now(tzone = "UTC"), lat = 0, lon = 0, eq.of.time)
#'
julian.day <- function(time) {
  2440587.79166667 + as.numeric(julian(time))
}

#' @rdname julian.day
#'
julian.century <- function(time) {
  (2440587.79166667 + as.numeric(julian(time)) - 2451545) / 36525
}

#' @rdname julian.day
#'
geom_mean_lon_sun <- function(x) {
  (280.46646 + x * (36000.76983 + x * 0.0003032)) %% 360
}

#' @rdname julian.day
#'
geom_mean_anom_sun <- function(x) {
  357.52911 + x * (35999.05029 - 0.0001537 * x)
}

#' @rdname julian.day
#'
eccent_earth_orbit <- function(x) {
  0.016708634 - x * (0.000042037 + 0.0000001267 * x)
}

#' @rdname julian.day
#'
sun_eq_of_ctr <- function(x, anom) {
  anom.rad <- anom / 180 * pi
  sin(anom.rad) * (1.914602 - x * (0.004817 + 0.000014 * x )) +
    sin(2 * anom.rad) * (0.019993 - 0.000101 * x) +
    sin(3 * anom.rad) * 0.000289
}

#' @rdname julian.day
#'
sun_rad_vector <- function(eccent, anom) {
  anom.rad <- anom / 180 * pi
  (1.000001018*  (1 - eccent^2))/(1 + eccent * cos(anom.rad))
}

#' @rdname julian.day
#'
sun_app_lon <- function(x, lon) {
  lon - 0.00569 - 0.00478 * sin((125.04 - 1934.136 * x) / 180 * pi)
}

#' @rdname julian.day
#'
mean_obliq_eclip <- function(x) {
  23 + (26 + ((21.448 -
                 x * (46.815 + x * (0.00059 - x * 0.001813)))) / 60) / 60
}

#' @rdname julian.day
#'
obliq_corr <- function(x, eclip) {
  eclip + 0.00256 * cos((125.04 - 1934.136 * x) / 180 * pi)
}

#' @rdname julian.day
#'
sun_rt_ascen <- function(app.lon, obliq.corr) {
  app.lon.rad <- app.lon / 180 * pi
  obliq.corr.rad <- obliq.corr / 180 * pi
  atan2(cos(obliq.corr.rad) * sin(app.lon.rad),
        cos(app.lon.rad )) / pi * 180
}

#' @rdname julian.day
#'
sun_declin <- function(app.lon, obliq.corr) {
  app.lon.rad <- app.lon / 180 * pi
  obliq.corr.rad <- obliq.corr / 180 * pi
  asin(sin(obliq.corr.rad) * sin(app.lon.rad)) /
    pi * 180
}

#' @rdname julian.day
#'
var_y <- function(obliq.corr) {
  x <- obliq.corr / 360 * pi
  tan(x)^2
}

#' @rdname julian.day
#'
eq_of_time <- function(mean.lon,
                       eccent.earth,
                       anom.mean,
                       var.y) {
  mean.lon.rad <- mean.lon / 180 * pi
  anom.mean.rad <- anom.mean / 180 * pi
  4 * ((var.y * sin(2 * mean.lon.rad) -
         2 * eccent.earth * sin(anom.mean.rad) +
         4 * eccent.earth * var.y * sin(anom.mean.rad) * cos(2 * mean.lon.rad) -
         0.5 * var.y^2 * sin(4 * mean.lon.rad) -
         1.25 * eccent.earth^2 * sin(2 * anom.mean.rad)) / pi * 180)
}

#' @rdname julian.day
#'
ha_sunrise <- function(lat, declin, ang = 0) {
  lat.rad <- lat / 180 * pi
  declin.rad <- declin / 180 * pi
  acos(cos((90 + ang) / 180 * pi) /            # 90.833 in NOAAs code
         (cos(lat.rad) * cos(declin.rad)) -
         tan(lat.rad) * tan(declin.rad))  / pi * 180
}

#' @rdname julian.day
#'
solar_noon <- function(lon, eq.of.time) {
  (720 - 4 * lon - eq.of.time) / 1440
}

#' @rdname julian.day
#'
sunrise <- function(noon, ha.sunrise) {
  (noon * 1440 -
     ha.sunrise * 4) / 1440
}

#' @rdname julian.day
#'
sunset <- function(noon, ha.sunrise) {
  (noon * 1440 +
     ha.sunrise * 4) / 1440
}

#' @rdname julian.day
#'
sunlight_duration <- function(ha.sunrise, unit.out = "hours") {
  8 * ha.sunrise
}

#' @rdname julian.day
#'
true_solar_time <- function(time, lat, lon, eq.of.time) {
 time + lubridate::seconds((eq.of.time + 4 * lon) * 60)
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
    when <- as_datetime(when, tz = tz.target)
  }
  (as.numeric(force_tz(when, tz.reference)) -
    as.numeric(force_tz(when, tz.target))) / 3600
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
#' @return \code{day_night} returns a data.fraame with variables sunrise time,
#'   sunset time, day length, night length. Each element of the list is a vector
#'   of the same length as the argument supplied for date.
#'
#' @note If twilight is a numeric vector of length two, the element with index 1
#'   is used for sunrise and that with index 2 for sunset.
#'
#' @family astronomy related functions
#'
#' @name day_night
#' @export
#' @examples
#' library(lubridate)
#' my.geocode <- data.frame(lat = 60, lon = 25)
#' day_night2(ymd("2015-05-30"), geocode = my.geocode, twilight = "civil")
#'
day_night <- function(date = lubridate::today(),
                      tz = Sys.timezone(),
                      geocode = data.frame(lon = 0, lat = 0),
                      twilight = "none",
                      unit.out = "hour") {
  stopifnot(! anyNA(date))
  stopifnot(is.data.frame(geocode))

  # if date is a vector or list for convenience we vectorize
  if (length(date) > 1) {
    first.iter <- TRUE
    for (d in date) {
      zz <- day_night(date = d,
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

  # from here onwards we are dealing with a single date

  if (unit.out %in% c("datetime", "date")) {
    duration.unit.out <- "hour"
  } else {
    duration.unit.out <- unit.out
  }

  twilight.angles <- twilight2angle(twilight)

  date <- lubridate::as_date(date, tz = tz)

  noon.of.date <- lubridate::as_datetime(date) + lubridate::hours(12)
  cent <- julian.century(noon.of.date)

  tz.diff <- tz_time_diff(noon.of.date, tz.target = tz)

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

  # up to here the calcualtions are not dependent on geographic location
  # now we vectorize for geocodes

  if (length(tz) == 1 && nrow(geocode) > 1) {
    tz <- rep(tz, nrow(geocode))
  }

  stopifnot(length(tz) == nrow(geocode))

  first.iter <- TRUE
  for (i in 1:nrow(geocode)) {
    lon <- geocode[i, "lon"]
    lat <- geocode[i, "lat"]

    # we do the calculations

    solar.noon <- solar_noon(lon, eq.of.time)

    if (twilight.angles[1] == twilight.angles[2]) {
      ha.sunrise <- ha_sunrise(lat, sun.declin, ang = -twilight.angles[1])
      sunrise <- sunrise(solar.noon, ha.sunrise)
      sunset <- sunset(solar.noon, ha.sunrise)
    } else {
      ha.sunrise1 <- ha_sunrise(lat, sun.declin, ang = -twilight.angles[1])
      ha.sunrise2 <- ha_sunrise(lat, sun.declin, ang = -twilight.angles[2])
      sunrise <- sunrise(solar.noon, ha.sunrise1)
      sunset <- sunset(solar.noon, ha.sunrise2)
    }

    sunrise.time <- lubridate::as_datetime(date, tz = tz[i]) +
      lubridate::seconds(sunrise * 86400)
    noon.time    <- lubridate::as_datetime(date, tz = tz[i]) +
      lubridate::seconds(solar.noon * 86400)
    sunset.time  <- lubridate::as_datetime(date, tz = tz[i]) +
      lubridate::seconds(sunset * 86400)

    daylength.duration <- sunset.time - sunrise.time
    daylength.hours <- as.numeric(daylength.duration)

    # we assemble the data frame for one geographic location

    if (unit.out %in% c("date", "datetime")) {
      yy <- data.frame(day           = date,
                       tz            = tz[i],
                       twilight.rise = twilight[1],
                       twilight.set  = ifelse(length(twilight) == 1,
                                            twilight[1], twilight[2]),
                       longitude     = lon,
                       latitude      = lat,
                       sunrise       = with_tz(sunrise.time, tzone = tz[i]),
                       noon          = with_tz(noon.time, tzone = tz[i]),
                       sunset        = with_tz(sunset.time, tzone = tz[i]),
                       daylength     = daylength.hours,
                       nightlength   = 24 - daylength.hours
      )
    } else if (unit.out %in% c("hours", "hour")) {
      yy <- data.frame(day           = date,
                       tz            = tz[i],
                       twilight.rise = twilight[1],
                       twilight.set  = ifelse(length(twilight) == 1,
                                                   twilight[1], twilight[2]),
                       longitude     = lon,
                       latitude      = lat,
                       sunrise       = (sunrise * 24 + tz.diff) %% 24,
                       noon          = (solar.noon * 24 + tz.diff) %% 24,
                       sunset        = (sunset * 24 + tz.diff) %% 24,
                       daylength     = daylength.hours,
                       nightlength   = 24 - daylength.hours
      )
    }

    # we bind the data frames together
    if (first.iter) {
      y <- yy
      first.iter <- FALSE
    } else {
      y <- rbind(y, yy)
    }
  }
  y
}



