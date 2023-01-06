#' Solar angles
#'
#' Function \code{sun_angles()} returns the solar angles and Sun to Earth
#' relative distance for given times and locations using a very precise
#' algorithm. Convenience functions \code{sun_azimuth()},
#' \code{sun_elevation()}, \code{sun_zenith_angle()} and
#' \code{distance_to_sun()} are wrappers on \code{sun_angles()} that return
#' individual vectors.
#'
#' @param time A "vector" of POSIXct Time, with any valid time zone (TZ) is
#'   allowed, default is current time.
#' @param tz character string indicating time zone to be used in output.
#' @param geocode data frame with variables lon and lat as numeric values
#'   (degrees), nrow > 1, allowed.
#' @param use.refraction logical Flag indicating whether to correct for
#'   fraction in the atmosphere.
#'
#' @return A \code{data.frame} with variables \code{time} (in same TZ as input),
#'   \code{TZ}, \code{solartime}, \code{longitude}, \code{latitude},
#'   \code{address}, \code{azimuth}, \code{elevation}, \code{declination},
#'   \code{eq.of.time}, \code{hour.angle}, and \code{distance}. If a data frame
#'   with multiple rows is passed to \code{geocode} and a vector of times longer
#'   than one is passed to \code{time}, sun position for all combinations of
#'   locations and times are returned by \code{sun_angles}. Angles are expressed
#'   in degrees, \code{solartime} is a vector of class \code{"solar.time"},
#'   \code{distance} is expressed in relative sun units.
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
#'   latitudes 90 and -90 degrees, i.e. exactly at the poles. The function is
#'   vectorized and in particular passing a vector of times for a single geocode
#'   enhances performance very much as the equation of time, the most time
#'   consuming step, is computed only once.
#'
#'   For improved performance, if more than one angle is needed it
#'   is preferable to directly call \code{sun_angles} instead of the wrapper
#'   functions as this avoids the unnecesary recalculation.
#'
#' @section Important!: Given an instant in time and a time zone, the date is
#'   computed from these, and may differ by one day to that at the location
#'   pointed by \code{geocode} at the same instant in time, unless the argument
#'   passed to \code{tz} matches the time zone at this location.
#'
#' @note There exists a different R implementation of the same algorithms called
#'   "AstroCalcPureR" available as function \code{astrocalc4r} in package
#'   'fishmethods'. Although the equations used are almost all the same, the
#'   function signatures and which values are returned differ. In particular,
#'   the present implementation splits the calculation into two separate
#'   functions, one returning angles at given instants in time, and a separate
#'   one returning the timing of events for given dates.
#'
#' @references
#' The primary source for the algorithm used is the book:
#' Meeus, J. (1998) Astronomical Algorithms, 2 ed., Willmann-Bell, Richmond,
#' VA, USA. ISBN 978-0943396613.
#'
#' A different implementation is available at
#' \url{https://github.com/NEFSC/READ-PDB-AstroCalc4R/}.
#'
#' An interactive web page using the same algorithms is available at
#' \url{https://gml.noaa.gov/grad/solcalc/}. There are small
#' differences in the returned times compared to our function that seem to be
#' related to the estimation of atmospheric refraction (about 0.1 degrees).
#'
#' @export
#'
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
  for (i in seq_len(nrow(geocode))) {
    temp <- sun_angles_fast(time = time,
                            tz = tz,
                            geocode = dplyr::slice(geocode, i),
                            use.refraction = use.refraction)
    z[[i]] <- temp # needed so that class attribute is retained
  }
  # we supress warning of dropped attributes and restore them
  z <- suppressWarnings(dplyr::bind_rows(z))
  class(z[["solartime"]]) <- class(temp[["solartime"]])

  # we could use rbind instead of dplyr::bind_rows as the second drops the class attribute
  # for solartime.
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
  sun.anom <- sun.anom.mean + delta
  sun.dist <- sun_rad_vector(eccent.earth, sun.anom)
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
  solar.time <- solar.time  %% 24 # needed for DST
  class(solar.time) <- c("solar_time", class(solar.time))

  tibble::tibble(time = lubridate::with_tz(time, tz),
                 tz = rep(tz, length(time)),
                 solartime = solar.time,
                 longitude = rep(lon, length(time)),
                 latitude = rep(lat, length(time)),
                 address = rep(address, length(time)),
                 azimuth = azimuth.angle,
                 elevation = elevation.angle,
                 declination = sun.declin,
                 eq.of.time = eq.of.time,
                 hour.angle = hour.angle,
                 distance = sun.dist,
                 .name_repair = "minimal")
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

#' @rdname sun_angles
#'
#' @export
#'
distance_to_sun <- function(time = lubridate::now(),
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
                  use.refraction = use.refraction)[["distance"]]
}

#' Time difference between two time zones
#'
#' Returns the difference in local time expressed in hours between two time
#' zones at a given instant in time. The difference due to daylight saving time
#' or Summer and Winter time as well as historical changes in time zones are
#' taken into account.
#'
#' @note This function is implemented using functions from package 'lubridate'.
#'   For details on the handling of time zones, please, consult the
#'   documentation for \code{\link{Sys.timezone}} about system differences in
#'   time zone names and handling.
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
#' Functions for calculating the timing of solar positions, given geographical
#' coordinates and dates. They can be also used to find the time for an
#' arbitrary solar elevation between 90 and -90 degrees by supplying "twilight"
#' angle(s) as argument.
#'
#' @param date "vector" of \code{POSIXct} times or\code{Date} objects, any valid
#'   TZ is allowed, default is current date at Greenwich matching the default
#'   for \code{geocode}.
#' @param tz character vector indicating time zone to be used in output and to
#'   interpret \code{Date} values passed as argument to \code{date}.
#' @param geocode data frame with one or more rows and variables lon and lat as
#'   numeric values (degrees). If present, address will be copied to the output.
#' @param twilight character string, one of "none", "rim", "refraction",
#'   "sunlight", "civil", "nautical", "astronomical", or a \code{numeric} vector
#'   of length one, or two, giving solar elevation angle(s) in degrees (negative
#'   if below the horizon).
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
#'   0 degrees. "rim": upper rim of solar disk at the horizon or solar elevation
#'   = -0.53 / 2. "refraction": solar elevation = 0 degrees + refraction
#'   correction. "sunlight": upper rim of solar disk corrected for refraction,
#'   which is close to the value used by the online NOAA Solar Calculator.
#'   "civil": -6 degrees, "naval": -12 degrees, and "astronomical": -18 degrees.
#'   Unit names for output are as follows: "day", "hours", "minutes" and
#'   "seconds" times for sunrise and sunset are returned as times-of-day since
#'   midnight expressed in the chosen unit. "date" or "datetime" return the same
#'   times as datetime objects with TZ set (this is much slower than "hours").
#'   Day length and night length are returned as numeric values expressed in
#'   hours when `"datetime"' is passed as argument to \code{unit.out}. If
#'   twilight is a numeric vector of length two, the element with index 1 is
#'   used for sunrise and that with index 2 for sunset.
#'
#'   \code{is_daytime()} supports twilight specifications by name, a test
#'   like \code{sun_elevation() > 0} may be used directly for a numeric angle.
#'
#' @seealso \code{\link{sun_angles}}.
#'
#' @note Function \code{day_night()} is an implementation of Meeus equations as
#'   used in NOAAs on-line web calculator, which are very precise and valid for
#'   a very broad range of dates. For sunrise and sunset the times are affected
#'   by refraction in the atmosphere, which does in turn depend on weather
#'   conditions. The effect of refraction on the apparent position of the sun is
#'   only an estimate based on "typical" conditions. The more tangential to the
#'   horizon is the path of the sun, the larger the effect of refraction is on
#'   the times of visual occlusion of the sun behind the horizon---i.e. the
#'   largest timing errors occur at high latitudes. The computation is not
#'   defined for latitudes 90 and -90 degrees, i.e. at the poles.
#'
#'   There exists a different R implementation of the same algorithms called
#'   "AstroCalcPureR" available as function \code{astrocalc4r} in package
#'   'fishmethods'. Although the equations used are almost all the same, the
#'   function signatures and which values are returned differ. In particular,
#'   the implementation in 'photobiology' splits the calculation into two
#'   separate functions, one returning angles at given instants in time, and a
#'   separate one returning the timing of events for given dates. In
#'   'fishmethods' (= 1.11-0) there is a bug in function astrocalc4r() that
#'   affects sunrise and sunset times. The times returned by the functions in
#'   package 'photobiology' have been validated against the NOAA base
#'   implementation.
#'
#'   In the current implementation functions \code{sunrise_time},
#'   \code{noon_time}, \code{sunset_time}, \code{day_length},
#'   \code{night_length} and \code{is_daytime} are all wrappers
#'   on \code{day_night}, so if more than one quantity is needed it is
#'   preferable to directly call \code{day_night} and extract the different
#'   components from the returned list.
#'
#' @section Warning: Be aware that R's \code{Date} class does not save time zone
#'   metadata. This can lead to ambiguities in the current implementation
#'   based on time instants. The argument passed to \code{date} should be
#'   of class \code{POSIXct}, in other words an instant in time, from which
#'   the correct date will be computed based on the \code{tz} argument.
#'
#'   The time zone in which times passed to \code{date} as argument are
#'   expressed does not need to be the local one or match the geocode, however,
#'   the returned values will be in the same time zone as the input.
#'
#' @references
#' The primary source for the algorithm used is the book:
#' Meeus, J. (1998) Astronomical Algorithms, 2 ed., Willmann-Bell, Richmond,
#' VA, USA. ISBN 978-0943396613.
#'
#' A different implementation is available at
#' \url{https://apps-nefsc.fisheries.noaa.gov/AstroCalc4R/} and in R paclage
#' 'fishmethods'. In 'fishmethods' (= 1.11-0) there is a bug in function
#' astrocalc4r() that affects sunrise and sunset times.
#'
#' An interactive web page using the same algorithms is available at
#' \url{https://gml.noaa.gov/grad/solcalc/}. There are small
#' differences in the returned times compared to our function that seem to be
#' related to the estimation of atmospheric refraction (about 0.1 degrees).
#'
#' @export
#' @examples
#' library(lubridate)
#'
#' my.geocode <- data.frame(lon = 24.93838,
#'                          lat = 60.16986,
#'                          address = "Helsinki, Finland")
#'
#' day_night(ymd("2015-05-30", tz = "EET"),
#'           geocode = my.geocode)
#' day_night(ymd("2015-05-30", tz = "EET") + days(1:10),
#'           geocode = my.geocode,
#'           twilight = "civil")
#' sunrise_time(ymd("2015-05-30", tz = "EET"),
#'              geocode = my.geocode)
#' noon_time(ymd("2015-05-30", tz = "EET"),
#'           geocode = my.geocode)
#' sunset_time(ymd("2015-05-30", tz = "EET"),
#'             geocode = my.geocode)
#' day_length(ymd("2015-05-30", tz = "EET"),
#'            geocode = my.geocode)
#' day_length(ymd("2015-05-30", tz = "EET"),
#'            geocode = my.geocode,
#'            unit.out = "day")
#' is_daytime(ymd("2015-05-30", tz = "EET") + hours(c(0, 6, 12, 18, 24)),
#'            geocode = my.geocode)
#' is_daytime(ymd_hms("2015-05-30 03:00:00", tz = "EET"),
#'            geocode = my.geocode)
#' is_daytime(ymd_hms("2015-05-30 00:00:00", tz = "UTC"),
#'            geocode = my.geocode)
#' is_daytime(ymd_hms("2015-05-30 03:00:00", tz = "EET"),
#'            geocode = my.geocode,
#'            twilight = "civil")
#' is_daytime(ymd_hms("2015-05-30 00:00:00", tz = "UTC"),
#'            geocode = my.geocode,
#'            twilight = "civil")
#'
day_night <- function(date = lubridate::now(tzone = "UTC"),
                      tz = ifelse(lubridate::is.Date(date),
                                  "UTC",
                                  lubridate::tz(date)),
                      geocode = tibble::tibble(lon = 0,
                                               lat = 51.5,
                                               address = "Greenwich"),
                      twilight = "none",
                      unit.out = "hours") {
  stopifnot(!anyNA(date)) # NAs could be propagated instead
  tz <- unique(tz)
  if (length(tz) > 1L) {
    tz <- tz[1]
    warning("'tz' is a heterogeneous vector, using only: ", tz)
  }
  geocode <- validate_geocode(geocode)
  if (any(lubridate::is.Date(date))) {
    date <- as.POSIXct(date, tz = tz)
  }
  # as 'date' is not a Date but a time, we find the corresponding date in UTC
  # as calculations are done in UTC time
  date <- lubridate::with_tz(date, tzone = "UTC")
  date <- lubridate::floor_date(date, unit = "days")

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

  z <- list()
  for (i in seq_len(nrow(geocode))) {
    temp <- day_night_fast(date = date,
                           tz = tz, # used for returned times
                           geocode = dplyr::slice(geocode, i),
                           twilight = twilight,
                           unit.out = unit.out)
    z[[i]] <- temp
  }

  z <- do.call(rbind, z)

  # we use rbind instead of dplyr::bind_rows as the second drops the class attribute
  # for solartime.

  # assertion
  if (any(!is.na(z[["daylength"]]) & z[["daylength"]] < 0) ||
      any(!is.na(z[["nightlength"]]) & z[["nightlength"]] < 0)) {
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
  # date should be always a POSIXct object with tz set to "UTC" at hms all set to zero.
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
                       day = 1/24,
                       1) # default

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

  noon.of.date <- lubridate::with_tz(date, tzone = "UTC") + 43200 # faster

  cent <- julian_century(noon.of.date)

  tz.diff <- tz_time_diff(noon.of.date, tz.target = tz)

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
    sunrise.time <- date +
      lubridate::seconds(sunrise * 86400)
    noon.time    <- date +
      lubridate::seconds(solar.noon * 86400)
    sunset.time  <- date +
      lubridate::seconds(sunset * 86400)

    tibble::tibble(day           = date,
                   tz            = rep(!!tz, length(date)),
                   twilight.rise = rep(twilight.angles[1], length(date)),
                   twilight.set  = rep(twilight.angles[2], length(date)),
                   longitude     = rep(lon, length(date)),
                   latitude      = rep(lat, length(date)),
                   address       = rep(address, length(date)),
                   sunrise       = lubridate::with_tz(sunrise.time, tzone = !!tz),
                   noon          = lubridate::with_tz(noon.time, tzone = !!tz),
                   sunset        = lubridate::with_tz(sunset.time, tzone = !!tz),
                   daylength     = daylength.hours,
                   nightlength   = 24 - daylength.hours,
                   .name_repair  = "minimal"
    )
  } else if (unit.out %in% c("day", "hour", "minute", "second")) {
    sunrise.tod <- (sunrise * 24 + tz.diff) %% 24
    noon.tod <- (solar.noon * 24 + tz.diff) %% 24
    sunset.tod <- (sunset * 24 + tz.diff) %% 24

    tibble::tibble(day           = date,
                   tz            = rep(!!tz, length(date)),
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
}

#' @rdname day_night
#'
#' @return \code{is_daytime()} returns a logical vector, with \code{TRUE} for
#'   day time and \code{FALSE} for night time.
#'
#' @export
#'
is_daytime <- function(date = lubridate::now(tzone = "UTC"),
                       tz = ifelse(lubridate::is.Date(date),
                                   "UTC",
                                   lubridate::tz(date)),
                       geocode = tibble::tibble(lon = 0,
                                                lat = 51.5,
                                                address = "Greenwich"),
                       twilight = "none",
                       unit.out = "hours") {
  if (!lubridate::is.POSIXct(date)) {
    warning("'date' must be a 'POSIXct' vector")
    return(rep(NA, length(date)))
  }
  z <- day_night(date = date,
                 tz = tz,
                 geocode = geocode,
                 twilight = twilight,
                 unit.out = "datetime")

  date > z[["sunrise"]] & date < z[["sunset"]]
}


#' twilight argument check and conversion
#'
#' @return numeric Solar elevation angle at sunrise or sunset
#' @keywords internal
#'
twilight2angle <- function(twilight) {
  # refraction depends on the elevation angle (ha), atmospheric pressure and temperature
  # we assume atmospheric pressure 1010 hPa and temperature 10 C
  # constants below are approximate
  # function atm_refraction_approx() can be used to estimate refraction for other angles
  if (length(twilight) == 0) {
    warning("'twilight' too short or NULL")
    twilight_angle <- c(NA_real_, NA_real_)
  } else if (is.character(twilight)) {
    if (length(twilight) != 1) {
      warning("'twilight' is character but of length > 1, using first value")
      twilight <- twilight[1L]
    }
    twilight_angle <-
      switch(twilight,
             none = c(0, 0), # center of solar disk
             rim = c(-0.53/2, -0.53/2), # upper rim of solar disk
             refraction = c(-0.4819444, -0.4819444), # center of solar disk, refraction corrected
             sunlight = c(-0.833, -0.833),
             civil = c(-6, -6),
             nautical = c(-12, -12),
             astronomical = c(-18, -18),
             {warning("Unrecognised 'twilight' value: ", twilight); c(NA_real_, NA_real_)}
      )
  } else if (is.numeric(twilight)) {
    if (length(twilight) == 1) {
      twilight_angle <- rep(twilight, 2)
    } else if (length(twilight) == 2) {
      twilight_angle <- twilight
    } else {
      stop("'twilight' vector longer than 2")
    }
    if (any(!is.na(twilight_angle) &
            twilight_angle < -90 | twilight_angle > 90 ) ) {
      warning("Off range twilight angle(s) replaced with NAs.")
      twilight_angle <-
        ifelse(!is.na(twilight_angle) &
                 twilight_angle > 90 | twilight_angle < -90,
               NA_real_, twilight_angle)
    }
  } else {
    warning("'twilight' must be numeric or character, but is ", class(twilight))
    twilight_angle <- c(NA_real_, NA_real_)
  }
  twilight_angle
}

#' @rdname day_night
#' @export
#' @return \code{noon_time}, \code{sunrise_time} and \code{sunset_time} return a
#'   vector of POSIXct times
#'
noon_time <- function(date = lubridate::now(tzone = "UTC"),
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
#'
sunrise_time <- function(date = lubridate::now(tzone = "UTC"),
                         tz = lubridate::tz(date),
                         geocode = tibble::tibble(lon = 0,
                                                  lat = 51.5,
                                                  address = "Greenwich"),
                         twilight = "sunlight",
                         unit.out = "datetime") {
 #  stopifnot(length(date) == 1L || nrow(geocode) == 1L)
  day_night(date = date,
            tz = tz,
            geocode = geocode,
            twilight = twilight,
            unit.out = unit.out)[["sunrise"]]
}

#' @rdname day_night
#'
#' @export
#'
sunset_time <- function(date = lubridate::now(tzone = "UTC"),
                        tz = lubridate::tz(date),
                        geocode = tibble::tibble(lon = 0,
                                                 lat = 51.5,
                                                 address = "Greenwich"),
                        twilight = "sunlight",
                        unit.out = "datetime") {
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
day_length <- function(date = lubridate::now(tzone = "UTC"),
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
night_length <- function(date = lubridate::now(tzone = "UTC"),
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

#' Convert datetime to time-of-day
#'
#' Convert a datetime into a time of day expressed in hours, minutes or seconds
#' from midnight in local time for a time zone. This conversion is useful when
#' time-series data for different days needs to be compared or plotted based on
#' the local time-of-day.
#'
#' @param x a datetime object accepted by lubridate functions
#' @param unit.out character string, One of "tod_time", "hours", "minutes", or "seconds".
#' @param tz character string indicating time zone to be used in output.
#'
#' @return A numeric vector of the same length as \code{x}. If
#'   \code{unit.out = "tod_time"} an object of class \code{"tod_time"} which
#'   the same as for \code{unit.out = "hours"} but with the class attribute
#'   set, which dispatches to special \code{format()} nad \code{print()}
#'   methods.
#'
#' @seealso \code{\link{solar_time}}
#'
#' @family Time of day functions
#'
#' @export
#'
#'
#' @examples
#' library(lubridate)
#' my_instants <- ymd_hms("2020-05-17 12:05:03") + days(c(0, 30))
#' my_instants
#' as_tod(my_instants)
#' as_tod(my_instants, unit.out = "tod_time")
#'
as_tod <- function(x, unit.out = "hours", tz = NULL) {
  stopifnot(lubridate::is.timepoint(x))
  if (!is.null(tz)) {
    x <- lubridate::with_tz(x, tzone = tz[1])
  }
  if (unit.out == "tod_time") {
    tod <- lubridate::hour(x) + lubridate::minute(x) / 60 + lubridate::second(x) / 3600
    class(tod) <- c("tod_time", class(tod))
    tod
  } else if (unit.out == "hours") {
    lubridate::hour(x) + lubridate::minute(x) / 60 + lubridate::second(x) / 3600
  } else if (unit.out == "minutes") {
    lubridate::hour(x) * 60 + lubridate::minute(x) + lubridate::second(x) / 60
  } else if (unit.out == "seconds") {
    lubridate::hour(x) * 3600 + lubridate::minute(x) * 60 + lubridate::second(x)
  } else {
    stop("Unrecognized 'unit.out': ", unit.out)
  }
}

#' Encode in a Common Format
#'
#' Format a \code{tod_time} object for pretty printing
#'
#' @param x an R object
#' @param ... ignored
#' @param sep character used as separator
#'
#' @family Time of day functions
#'
#' @export
#'
format.tod_time <- function(x, ..., sep = ":") {
  hours <- as.integer(trunc(x))
  minutes <- as.integer((x * 60) %% 60)
  seconds <- as.integer((x * 3600) %% 60)
  fmt <- paste(rep("%02d", 3), collapse = sep)
  time_string <-
    sprintf(fmt = fmt, hours, minutes, seconds)
  time_string
}

#' Print time-of-day objects
#'
#' @param x an R object
#' @param ... passed to \code{format} method
#'
#' @family Time of day functions
#'
#' @note Default is to print the underlying \code{numeric} vector as a solar time.
#'
#' @export
#'
print.tod_time <- function(x, ...) {
  print(format(x, ...))
  invisible(x)
}

#' Local solar time
#'
#' \code{solar_time()} computes the time of day expressed in seconds since the
#' astronomical midnight using and instant in time and a geocode as input. Solar
#' time is useful when we want to plot data according to the local solar time
#' rather than the local time in use at a time zone. How the returned instant in
#' time is expressed depends on the argument passed to \code{unit.out}.
#'
#' @details Solar time is determined by the position of the sun in the sky and
#' it almost always differs from the time expressed in the local time
#' coordinates in use. The differences can vary from a few minutes up to a
#' couple of hours depending on the exact location within the time zone and the
#' use or not of daylight saving time.
#'
#' @param time POSIXct Time, any valid time zone (TZ) is allowed, default is
#'   current time
#' @param geocode data frame with variables lon and lat as numeric values
#'   (degrees).
#' @param unit.out character string, One of "datetime", "time", "hour", "minute", or
#'   "second".
#'
#' @seealso \code{\link{as_tod}}
#'
#' @family Local solar time functions
#'
#' @note The algorithm is approximate, it calculates the difference between
#'   local solar noon and noon in the time zone of \code{time} and uses this
#'   value for the whole day when converting times into solar time. Days are not
#'   exactly 24 h long. Between successive days the shift is only a few seconds,
#'   and this leads to a small jump at midnight.
#'
#' @section Warning!: Returned values are computed based on the time zone of the
#'   argument for parameter time. In the case of solar time, this timezone does
#'   not affect the result. However, in the case of solar dates the date part
#'   may be off by one day, if the time zone does not match the coordinates of
#'   the geocode value provided as argument.
#'
#' @return In all cases solar time is expressed as time since local astronomical
#'   midnight and, thus, lacks date information. If \code{unit.out = "time"}, a
#'   numeric value in seconds with an additional class attribute
#'   "solar_time"; if \code{unit.out = "datetime"}, a "POSIXct" value in seconds
#'   from midnight but with an additional class attribute "solar_date"; if
#'   \code{unit.out = "hour"} or \code{unit.out = "minute"} or \code{unit.out =
#'   "second"}, a numeric value.
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
#' @family Local solar time functions
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
#' @family Local solar time functions
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
#' @family Local solar time functions
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

#' Validate a geocode
#'
#' Test validity of a geocode or ensure that a geocode is valid.
#'
#' @details
#' \code{validate_geocode} Converts to tibble, checks data bounds, converts
#' address to character if it is not already a character vector, or add
#' character NAs if the address column is missing.
#'
#' \code{is_valid_geocode} Checks if a geocode is valid, returning 0L if not,
#' and the number of row otherwise.
#'
#' @param geocode data.frame with geocode data in columns \code{"lat"},
#'   \code{"lon"}, and possibly also \code{"address"}.
#'
#' @return A valid geocode stored in a tibble.
#'
#' @examples
#'
#' validate_geocode(NA)
#' validate_geocode(data.frame(lon = -25, lat = 66))
#'
#' is_valid_geocode(NA)
#' is_valid_geocode(1L)
#' is_valid_geocode(data.frame(lon = -25, lat = 66))
#'
#' na_geocode()
#'
#' @export
#'
validate_geocode <- function(geocode) {
  if (is.atomic(geocode) && (length(geocode) == 1L) && is.na(geocode)) {
    geocode <- na_geocode()
  } else if (is.data.frame(geocode)) {
    geocode <- tibble::as_tibble(geocode, .name_repair = "minimal")
  } else {
    stop("Bad geocode: ", format(geocode))
  }
  stopifnot(nrow(geocode) >= 1) # needs to be replace by generation of no output in all functions
  stopifnot(exists("lon", geocode), exists("lat", geocode))
  geocode[["lon"]] <- as.numeric(geocode[["lon"]]) # convert logical NA
  geocode[["lat"]] <- as.numeric(geocode[["lat"]]) # convert logical NA
  if (any(na.omit(geocode[["lon"]]) > 180 | na.omit(geocode[["lon"]]) < -180)) {
    stop("Longitude is off-range.")
  }
  if (any(na.omit(geocode[["lat"]]) > 89.99 | na.omit(geocode[["lat"]]) < -89.99)) {
    stop("Latitude is off-range.")
  }
  if (!exists("address", geocode)) {
    geocode[["address"]] <- NA_character_
  } else if (!is.character(geocode[["address"]])) {
    geocode[["address"]] <- as.character(geocode[["address"]])
  }
  geocode
}

#' @rdname validate_geocode
#'
#' @return FALSE for invalid, TRUE for valid.
#'
#' @export
#'
is_valid_geocode <- function(geocode) {
  if (!is.list(geocode)) return(FALSE)
  if (!is.data.frame(geocode)) {
    # walk list of geocodes using recursion
    is_valid <- all(sapply(geocode, is_valid_geocode))
  } else {
    is_valid <- nrow(geocode) >= 1L &&
      all(c("lon", "lat") %in% names(geocode)) &&
      all(c(is.numeric(geocode[["lon"]]), is.numeric(geocode[["lat"]]))) &&
      if ("address" %in% names(geocode)) is.character(geocode[["address"]]) else TRUE

    if (!is_valid && "address" %in% names(geocode) && is.factor(geocode[["address"]])) {
      warning("'address' is a factor instead of a character vector.")
    }
  }
  is_valid
}

#' @rdname validate_geocode
#'
#' @return FALSE for invalid, number of rows for valid.
#'
#' @export
#'
length_geocode <- function(geocode) {
  if (!is_valid_geocode(geocode)) return(NA_integer_)
  nrow(geocode)
}

#' @rdname validate_geocode
#'
#' @return A geo_code tibble with all fields set to suitable NAs.
#'
#' @export
#'
na_geocode <- function() {
  tibble::tibble(lon = NA_real_,
                 lat = NA_real_,
                 address = NA_character_)
}

