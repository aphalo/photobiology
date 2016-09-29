#' Solar angles
#'
#' This function returns the solar angles for a given time and location.
#'
#' @param time POSIXct Time, any valid time zone (TZ) is allowed, default is
#'   current time
#' @param geocode data frame with variables lon and lat as numeric values
#'   (degrees).
#' @param lon numeric Vector of longitudes (degrees) W is < 0, and E is > 0
#' @param lat numeric Vector of latitudes (degrees)
#' @param use_refraction logical Flag indicating whether to correct for
#'   fraction in the atmosphere
#'
#' @return A list with components time in same TZ as input, azimuth, elevation,
#'   diameter, and distance.
#'
#' @family astronomy related functions
#'
#' @references
#' Michalsky, J. J., 1988. "The Astronomical Almanac's algorithm for approximate
#' solar position (1950--2050)". Solar Energy, 227--235.
#'
#' Spencer, J. W., 1989. "Comments on The Astronomical Almanac's algorithm for
#' approximate solar position (1950--2050)." Solar Energy, 42, 353.
#'
#' Vignola, F.; Michalsky, J. & Stoffel, T., 2012 (Eds.) "Solar and infrared
#' radiation measurements." Boca Raton, CRC Press, ISBN 9781439851890.
#'
#' @note Several implementations of this algorithm are available in the
#'   internet. I have found FORTRAN, Perl and R versions. Not all these
#'   implementations correctly handle year 2000, which is an exception to the
#'   normal leap-year rule. They also differ on the handling of refraction. The
#'   algorithm used only accepts dates for years 1950 to 2050. A listing is also
#'   available in an appendix in Vignola et al. (2012). This implementation is
#'   not a direct copy or translation of any of these examples. The current
#'   version of the code owes much to Josh O'Brien's asnwer to a question in
#'   StackOverflow
#'   \url{http://stackoverflow.com/questions/8708048/position-of-the-sun-given-time-of-day-latitude-and-longitude}
#'
#' @export
#' @examples
#' require(lubridate)
#' sun_angles()
#' sun_angles(ymd_hms("2014-09-23 12:00:00"))
#' sun_angles(ymd_hms("2014-09-23 12:00:00"), lat=60, lon=0)
#'
sun_angles <- function(time = lubridate::now(),
                       geocode = NULL,
                       lon = 0, lat = 0,
                       use_refraction = FALSE)
{
  # validate arguments
  stopifnot(lubridate::is.POSIXct(time))
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  stopifnot(abs(lat) <= 90 + 1e-20)
  stopifnot(abs(lon) <= 180 + 1e-20)
  # take care of time zone
  tz <- lubridate::tz(time)
  t <- lubridate::with_tz(time, "UTC")
  # input can be a vector of times
  nt <- length(t)

  # if geocode argument supplied override lat and lon
  # locations can be also vectors
  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }
  # do recycling of arguments if needed
  nlon <- length(lon)
  nlat <- length(lat)
  # lat and lon of different lengths
  if (nlon == 1 && nlat > 1) {
    lon <- rep(lon, nlat)
    nlon <- nlat
  } else if (nlat == 1 && nlon > 1) {
    lat <- rep(lat, nlon)
    nlat <- nlon
  } else if (nlon != nlat) {
    stop("lengths of longitude and latitude must match")
  }
  # number of locations different from number of times
  if (nlon == 1 && nt > 1) {
    lon <- rep(lon, nt)
    nlon <- nt
    lat <- rep(lat, nt)
    nlat <- nt
  } else if (nt == 1 && nlon > 1) {
    t <- rep(t, nlon)
    nt <- nlon
  } else if (nt != nlon) {
    stop("lengths of t, latitude and longitude must match, or have length 1")
  }
  year <- lubridate::year(t)
  if (any(year < 1950) || any(year > 2050))
    stop("year=", year, " is outside accepted range")
  # this already corrects for leap years
  hour <- lubridate::hour(t) + lubridate::minute(t) / 60 + lubridate::second(t) / 3600
  time <- as.numeric(lubridate::as.duration(
    t - lubridate::ymd_hms("2000-01-01 12:00:00", tz = "UTC"))) / 3600 / 24
  # Ecliptic coordinates
  # Mean longitude
  mnlong <- 280.46 + 0.9856474 * time
  mnlong <- mnlong %% 360
  mnlong <- ifelse(mnlong < 0, mnlong + 360, mnlong)
  # Mean anomaly
  mnanom <- 357.528 + 0.9856003 * time
  mnanom <- mnanom %% 360
  mnanom <- ifelse(mnanom < 0, mnanom + 360, mnanom)
  rpd <- pi/180
  mnanom <- mnanom * rpd
  # Ecliptic longitude and obliquity of ecliptic
  eclong <- mnlong + 1.915 * sin(mnanom) + 0.02 * sin(2 * mnanom)
  eclong <- eclong %% 360
  eclong <- ifelse(eclong < 0, eclong + 360, eclong)
  oblqec <- 23.439 - 4e-07 * time
  eclong <- eclong * rpd
  oblqec <- oblqec * rpd
  # Celestial coordinates
  # Right ascension and declination
  num <- cos(oblqec) * sin(eclong)
  den <- cos(eclong)
  ra <- atan(num / den)
  ra <- ifelse(den < 0, ra + pi, ifelse(num < 0, ra + 2 * pi, ra))
  dec <- asin(sin(oblqec) * sin(eclong))
  # Local coordinates
  # Greenwich mean sidereal time
  # $h = $hour + $min / 60 + $sec / 3600;
  gmst <- 6.697375 + 0.0657098242 * time + hour
  gmst <- gmst %% 24
  gmst <- ifelse(gmst < 0, gmst + 24, gmst)
  # Local mean sidereal time
  lmst <- gmst + lon / 15
  lmst <- lmst %% 24
  lmst <- ifelse(lmst < 0, lmst + 24, lmst)
  lmst <- lmst * 15 * rpd
  # Hour angle
  ha <- lmst - ra
  ha <- ifelse(ha < (-pi), ha + 2 * pi, ha)
  ha <- ifelse(ha > pi, ha - 2 * pi, ha)
  # Latitude to radians
  lat <- lat * rpd
  # Solar zenith angle
  za <- acos(sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(ha))
  # Solar azimuth
  az <- acos(((sin(lat) * cos(za)) - sin(dec)) / (cos(lat) * sin(za)))
  # Solar elevation
  el <- asin(sin(dec) * sin(lat) + cos(dec) * cos(lat) * cos(ha))
  # Latitude to radians
  el <- el/rpd
  az <- az/rpd
  lat <- lat/rpd

  az <- ifelse(ha > 0, az + 180, 540 - az)
  az <- az %% 360
  # refraction correction
  if (use_refraction) {
    refrac <-
      ifelse(el >= 19.225, 0.00452 * 3.51823/tan(el * rpd),
             ifelse(el > (-0.766) & el < 19.225,
                    3.51823 * (0.1594 + el * (0.0196 + 2e-05 * el)) /
                      (1 + el * (0.505 + 0.0845 * el)),
                    0))
    el <- el + refrac
  }
  # solar distance and diameter
  soldst <- 1.00014 - 0.01671 * cos(mnanom) - 0.00014 * cos(2 * mnanom)
  soldia <- 0.5332 / soldst
  # assertion
  if (any(el < (-90)) || any(el > 90))
    stop("output el out of range")
  if (any(az < 0) || any(az > 360))
    stop("output az out of range")
  # return values
  return(list(time = lubridate::with_tz(t, tz),
              longitude = lon,
              latitude = lat,
              azimuth = az,
              elevation = el,
              diameter = soldia,
              distance = soldst))
}

#' Times for sun positions
#'
#' Functions for calculating the timing of solar positions by means of function
#' \code{sun_angles}, given geographical coordinates and dates. They can be also
#' used to find the time for an arbitrary solar elevation between 90 and -90
#' degrees by supplying "twilight" angle(s) as argument.
#'
#' @param date array of POSIXct times or Date objects, any valid TZ is allowed,
#'   default is current date
#' @param tz character string incading time zone to be used in output, default
#'   is system time zone
#' @param geocode data frame with variables lon and lat as numeric values
#'   (degrees).
#' @param lon numeric array of longitudes (degrees)
#' @param lat numeric array of latitudes (degrees)
#' @param twilight character string, one of "none", "civil", "nautical",
#'   "astronomical", or a \code{numeric} vector of length one, or two, giving
#'   solar elevation angle(s) in degrees (negative if below the horizon).
#' @param unit.out charater string, One of "date", "hour", "minute", or "second".
#'
#' @return \code{day_night} returns a list with fields sunrise time, sunset
#'   time, day length, night length. Each element of the list is a vector of the
#'   same length as the argument supplied for date.
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
#' day_length()
#' day_length(ymd("2015-05-30"), lat = 60, lon = 25)
#' day_length(ymd("2014-12-30"), lat = 60, lon = 25)
#' day_length(ymd("2015-05-30"), lat = 60, lon = 25, twilight = "civil")
#' sunrise_time(ymd("2015-05-30"), lat = 60, lon = 25, tz = "EET")
#' day_night(ymd("2015-05-30"), lat = 60, lon = 25, twilight = "civil")
#'
day_night <- function(date = lubridate::today(),
                      tz = "UTC",
                      geocode = NULL, lon = 0, lat = 0,
                      twilight = "none",
                      unit.out = "date") {
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }

  list(day         = as.Date(date),
       sunrise     = sunrise_time(date = date, tz = tz,
                                  lon = lon, lat = lat,
                                  twilight = twilight,
                                  unit.out = unit.out),
       noon        = noon_time(date = date, tz = tz,
                               lon = lon, lat = lat,
                               unit.out = unit.out),
       sunset      = sunset_time(date = date, tz = tz,
                                 lon = lon, lat = lat,
                                 twilight = twilight,
                                 unit.out = unit.out),
       daylength   = day_length(date = date, tz = tz,
                                lon = lon, lat = lat,
                                twilight = twilight,
                                unit.out = unit.out),
       nightlength = night_length(date = date, tz = tz,
                                  lon = lon, lat = lat,
                                  twilight = twilight,
                                  unit.out = unit.out) )
}

#' twilight argument check and conversion
#'
#' @return numeric  Solar elevation angle at sunrise or sunset
#' @keywords internal
twilight2angle <- function(twilight) {
  if (!is.numeric(twilight)) {
    if (twilight == "none") {
      twilight_angle <- c(0, 0)
    } else if (twilight == "civil") {
      twilight_angle <- c(-6, -6)
    } else if (twilight == "nautical") {
      twilight_angle <- c(-12, -12)
    } else if (twilight == "astronomical") {
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

#' date argument check and conversion
#'
#' @return numeric representtaion of the date
#' @keywords internal
date2seconds <- function(t, tz) {
  if (!lubridate::is.POSIXct(t)) {
    if (lubridate::is.instant(t)) {
      t <- as.POSIXct(t, tz = "UTC")
    } else {
      warning("t is not a valid time or date")
    }
  }
  t <- as.POSIXct(t, tz = tz)
  lubridate::hour(t) <- 0
  lubridate::minute(t) <- 0
  lubridate::second(t) <- 0
  as.numeric(t, tz = tz)
}

#' time argument check and conversion
#'
#' @return numeric representtaion of the date
#' @keywords internal
time2seconds <- function(t, tz) {
  if (!lubridate::is.POSIXct(t)) {
    if (lubridate::is.instant(t)) {
      t <- as.POSIXct(t, tz = "UTC")
    } else {
      warning("t is not a valid time or date")
    }
  }
  t <- as.POSIXct(t, tz = tz)
  as.numeric(t, tz = tz)
}

#' function to be numerically minimized
#'
#' @return an elevation angle delta
#' @keywords internal
altitude <- function(x, lon, lat, twlght_angl = 0){
  t_temp <- as.POSIXct(x, origin = lubridate::origin, tz = "UTC")
  return(sun_angles(t_temp,
                    lon = lon,
                    lat = lat)$elevation - twlght_angl)
}

#' @rdname day_night
#' @export
#' @return \code{noon_time}, \code{sunrise_time} and \code{sunset_time} return a
#'   vector of POSIXct times
noon_time <- function(date = lubridate::today(), tz = "UTC",
                      geocode = NULL, lon = 0, lat = 0,
                      twilight = NA, unit.out = "date") {
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }
  # we optimize for the simplest cases when vector arguments are all equal
  # which occurs frequently when noon time is called by solar_time()
  if (all(length(date) > 1, length(lon) > 1, length(lat) > 1) &&
      all.equal(length(date), length(lon), length(lat)) &&
      all.equal(date) &&
      all.equal(lon) &&
      all.equal(lat)) {
    date <- unique(date)
    lon <- unique(lon)
    lat <- unique(lat)
    length.out <- rep(1,length(date))
  } else if (length(date) > 1 && length(lon) == 1 && length(lat) == 1) {
    if (length(unique(date)) == 1) {
      length.out <- length(date)
      date <- unique(date)
    } else if (length(unique(date)) < length(date) && !is.unsorted(date)) {
      transitions <- c(which(diff(date) != 0), length(date))
      length.out <- c(transitions[1], diff(transitions))
      date <- unique(date)
    } else {
      length.out <- 1
    }
  } else {
    length.out <- 1
  }
  date_num <- sapply(date, date2seconds, tz = tz)
  times <- numeric()
  twlght_angl <- 0
  for (t_num in date_num) {
    noon <- try(
      stats::optimize(f = altitude,
                      interval = c(t_num + 7200, t_num + 86400 - 7200),
                      lon = lon, lat = lat,
                      twlght_angl = twlght_angl,
                      maximum = TRUE)$maximum
    )
    if (inherits(noon, "try-error")) {
      noon <- NA
    }
    times <- c(times, noon)
  }
  times <- as.POSIXct(times, tz = tz, origin = lubridate::origin)
  if (unit.out != "date") {
    times <- sapply(times, date2tod, unit.out = unit.out)
  }
  rep(times, length.out)
}

#' @rdname day_night
#'
#' @export
sunrise_time <- function(date = lubridate::today(), tz = "UTC",
                         geocode = NULL, lon = 0, lat = 0,
                         twilight = "none", unit.out = "date") {
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }
  noon <- noon_time(date = date, tz = tz,
                    lon = lon, lat = lat)
  noon_num <- sapply(noon, time2seconds, tz = tz)
  times <- numeric()
  twlght_angl <- twilight2angle(twilight)[1]
  for (t_num in noon_num) {
     rise <- try(
      stats::uniroot(f = altitude,
              lon = lon, lat = lat, twlght_angl = twlght_angl,
              lower = t_num - 86400/2, upper = t_num)$root,
      silent = TRUE)
    if (inherits(rise, "try-error")) {
      rise <- NA # never
    }
    times <- c(times, rise)
  }
  times <- as.POSIXct(times, tz = tz, origin = lubridate::origin)
  if (unit.out != "date") {
    times <- sapply(times, date2tod, unit.out = unit.out)
  }
  times
}

#' @rdname day_night
#' @export
#'
#' @note \code{night_length} returns the length of night-time conditions in one
#'   day (00:00:00 to 23:59:59), rather than the length of the night between two
#'   consequtive days.
sunset_time <- function(date = lubridate::today(), tz = "UTC",
                        geocode = NULL, lon = 0, lat = 0,
                        twilight = "none", unit.out = "date") {
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }
  noon <- noon_time(date = date, tz = tz,
                    lon = lon, lat = lat)
  noon_num <- sapply(noon, time2seconds, tz = tz)
  times <- numeric()
  twlght_angl <- twilight2angle(twilight)[2]
  for (t_num in noon_num) {
    set <- try(
      stats::uniroot(f = altitude,
              lon = lon, lat = lat,
              twlght_angl = twlght_angl,
              lower = t_num - 1, upper = t_num + 86400 / 2)$root,
      silent = TRUE)
    if (inherits(set, "try-error")) {
      set <- NA # never
    }
    times <- c(times, set)
  }
  times <- as.POSIXct(times, tz = tz, origin = lubridate::origin)
  if (unit.out != "date") {
    times <- sapply(times, date2tod, unit.out = unit.out)
  }
  times
}

#' @rdname day_night
#'
#' @export
#' @return \code{day_length} and \code{night_length} return numeric a vector
#'   giving the length in hours
day_length <- function(date = lubridate::today(), tz = "UTC",
                       geocode = NULL, lon = 0, lat = 0,
                       twilight = "none", unit.out = "hour") {
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }
  noon <- noon_time(date = date, tz = tz,
                    lon = lon, lat = lat)
  rise_time <- sunrise_time(date = date, tz = tz,
                            lon = lon, lat = lat,
                            twilight = twilight, unit.out = "date")
  set_time <- sunset_time(date = date, tz = tz,
                          lon = lon, lat = lat,
                          twilight = twilight, unit.out = "date")
  hours <- ifelse(is.na(rise_time) | is.na(set_time),
         ifelse(altitude(noon, lon = lon, lat = lat) > 0, 24, 0),
         set_time - rise_time)
  if (any(hours > 24) || any(hours < 0)) {
    hours <- NA_real_
    warning("Convergence failed!")
  }
  switch(unit.out,
         "date" = hours,
         "hour" = hours,
         "minute" = hours * 60,
         "second" = hours * 3600)
}

#' @rdname day_night
#'
#' @export
night_length <- function(date = lubridate::today(), tz = "UTC",
                         geocode = NULL, lon = 0, lat = 0,
                         twilight = "none", unit.out = "hour") {
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }
  hours <- 24 - day_length(date = date, tz = tz,
                           lon = lon, lat = lat,
                           twilight = twilight, unit.out = "hour")
  switch(unit.out,
         "date" = hours,
         "hour" = hours,
         "minute" = hours * 60,
         "second" = hours * 3600)
}

#' Convert date to time-of-day in hours
#'
#' @param date a date object accepted by lubridate functions
#'
#' @keywords internal
#'
date2tod <- function(date, unit.out) {
  if (unit.out == "hour") {
    lubridate::hour(date) + lubridate::minute(date) / 60 + lubridate::second(date) / 3600
  } else if (unit.out == "minute") {
    lubridate::hour(date) * 60 + lubridate::minute(date) + lubridate::second(date) / 60
  } else if (unit.out == "second") {
    lubridate::hour(date) * 3600 + lubridate::minute(date) * 60 + lubridate::second(date)
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
#' sol_d <- solar_date(lubridate::dmy_hms("21/06/2016 10:00:00", tz = "UTC"),
#'                     BA.geocode)
#' sol_d
#' class(sol_d)
#'
solar_time <- function(time = lubridate::now(), geocode = NULL)
{
  if (lubridate::is.POSIXlt(time)) {
    new <- as.POSIXct(time)
  } else {
    new <- time
  }
  solar.noon <- noon_time(date = lubridate::as_date(new),
                          tz = lubridate::tz(time),
                          geocode = geocode,
                          twilight = NA,
                          unit.out = "date")
  # solar time in hours from midnight
  solar.time <- as.double(new - (solar.noon - lubridate::hours(12)), units = "hours")
  attr(solar.time, "tzone") <- lubridate::tz(time)
  class(solar.time) <- c("solar_time", class(solar.time))
  solar.time
}

#' @rdname solar_time
#'
#' @return For \code{solar_date()} a date-time object with the class attr set
#'  to "solar.time". This is needed only for unambiguous formating and printing.
#'
#' @export
#'
solar_date <- function(time = lubridate::now(), geocode = NULL)
{
  solar.time <- solar_time(time = time, geocode = geocode)
  solar.date <-
    lubridate::floor_date(time, unit = "day") +
    lubridate::seconds(solar.time * 3600)
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
is.solar_time <- function(x) {
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
  tz <- attr(x, "tzone")
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
