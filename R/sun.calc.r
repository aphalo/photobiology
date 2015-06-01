#' Calculate solar angles
#'
#' This function returns the solar angles for a given time and location.
#'
#' @param t POSIXct Time, any valid time zone (TZ) is allowed, default is
#'   current time
#' @param lon numeric Vector of longitudes (degrees)
#' @param lat numeric Vector of latitudes (degrees)
#' @param use_refraction logical Flag indicating whether to correct for
#'   fraction in the atmosphere
#'
#' @return A list with components time in same TZ as input, azimuth, elevation,
#'   diameter, and distance.
#'
#' @keywords manip misc
#' @export
#' @examples
#' require(lubridate)
#' sun_angles()
#' sun_angles(ymd_hms("2014-09-23 12:00:00"))
#' sun_angles(ymd_hms("2014-09-23 12:00:00"), lat=60, lon=0)
#'
sun_angles <- function(t = now(), lon = 0, lat = 0, use_refraction = FALSE)
{
  if (!is.POSIXct(t)) {
    warning("Argument t is not a POSIXct time.")
    return(NA)
  }
  tz <- tz(t)
  t <- with_tz(t, "UTC")
  nt <- length(t)
  nlon <- length(lon)
  nlat <- length(lat)
  if (nlon != nlat)
    stop("lengths of longitude and latitude must match")
  if (nlon == 1) {
    lon <- rep(lon, nt)
    lat <- rep(lat, nt)
  }
  else {
    if (nt != nlon)
      stop("lengths of t, latitude and longitude must match, unless last two are of length 1")
  }
  year <- year(t)
  if (any(year < 1950) || any(year > 2050))
    stop("year=", year, " is outside acceptable range")
  day <- yday(t)
  if (any(day < 1) || any(day > 366))
    stop("day is not in range 1 to 366")
  hour <- hour(t) + minute(t) / 60 + second(t) / 3600
  if (any(hour < -13) || any(hour > 36))
    stop("hour outside range -13 to 36")
  if (any(lat < -90)) {
    warning("latitude(s) trimmed to range -90 to 90")
    lat[lat < -90] <- -90
  }
  if (any(lat > 90)) {
    warning("latitude(s) trimmed to range -90 to 90")
    lat[lat > 90] <- 90
  }
  if (any(lon < -180)) {
    warning("longitude(s) trimmed to range -180 to 180")
    lon[lon < -180] <- -180
  }
  if (any(lon > 180)) {
    warning("longitude(s) trimmed to range -180 to 180")
    lon[lon > 180] <- 180
  }
  delta <- year - 1949
  leap <- delta%/%4
  jd <- 32916.5 + (delta * 365 + leap + day) + hour / 24
  jd <- jd + ifelse(0 == (year%%100) & 0 != (year%%400), 1,
                    0)
  time <- jd - 51545
  mnlong <- 280.46 + 0.9856474 * time
  mnlong <- mnlong%%360
  mnlong <- mnlong + ifelse(mnlong < 0, 360, 0)
  mnanom <- 357.528 + 0.9856003 * time
  mnanom <- mnanom%%360
  mnanom <- mnanom + ifelse(mnanom < 0, 360, 0)
  rpd <- pi/180
  mnanom <- mnanom * rpd
  eclong <- mnlong + 1.915 * sin(mnanom) + 0.02 * sin(2 * mnanom)
  eclong <- eclong%%360
  eclong <- eclong + ifelse(eclong < 0, 360, 0)
  oblqec <- 23.439 - 4e-07 * time
  eclong <- eclong * rpd
  oblqec <- oblqec * rpd
  num <- cos(oblqec) * sin(eclong)
  den <- cos(eclong)
  ra <- atan(num/den)
  ra <- ra + ifelse(den < 0, pi, ifelse(num < 0, 2 * pi, 0))
  dec <- asin(sin(oblqec) * sin(eclong))
  gmst <- 6.697375 + 0.0657098242 * time + hour
  gmst <- gmst%%24
  gmst <- gmst + ifelse(gmst < 0, 24, 0)
  lmst <- gmst + lon/15
  lmst <- lmst%%24
  lmst <- lmst + ifelse(lmst < 0, 24, 0)
  lmst <- lmst * 15 * rpd
  ha <- lmst - ra
  ha <- ha + ifelse(ha < (-pi), 2 * pi, 0)
  ha <- ha - ifelse(ha > pi, 2 * pi, 0)
  el <- asin(sin(dec) * sin(lat * rpd) + cos(dec) * cos(lat *
                                                          rpd) * cos(ha))
  az <- asin(-cos(dec) * sin(ha)/cos(el))
  az <- ifelse(sin(dec) - sin(el) * sin(lat * rpd) > 0,
               ifelse(sin(az) < 0, az + 2 * pi, az), pi - az)
  el <- el/rpd
  az <- az/rpd
  if (use_refraction) {
    refrac <- ifelse(el >= 19.225, 0.00452 * 3.51823/tan(el * rpd),
                     ifelse(el > (-0.766) & el < 19.225,
                            3.51823 * (0.1594 + el * (0.0196 + 2e-05 * el))/(1 + el * (0.505 + 0.0845 * el)),
                            0))
    el <- el + refrac
  }
  soldst <- 1.00014 - 0.01671 * cos(mnanom) - 0.00014 * cos(2 * mnanom)
  soldia <- 0.5332/soldst
  if (any(el < (-90)) || any(el > 90))
    stop("output el out of range")
  if (any(az < 0) || any(az > 360))
    stop("output az out of range")
  return(list(time = with_tz(t, tz),
              azimuth = az,
              elevation = el,
              diameter = soldia,
              distance = soldst))
}

#' Calculate time of sunrise and sunset
#'
#' This function returns the times of sunrise and sunset for a given location
#' and date. It can be also used to find the time for an arbitrary solar
#' elevation between zenith and -30 degrees by supplying  a "twilight" angle
#' as argument.
#'
#' @param t array of POSIXct times, any valid TZ is allowed, default is current
#'   date
#' @param lon numeric array of longitudes (degrees)
#' @param lat numeric array of latitudes (degrees)
#' @param twilight character string, one of "none", "civil", "nautical",
#'   "astronomical", or a \code{numeric} value giving a solar elevation angle in
#'   degrees (negative if below the horizon)
#' @param tz character string incading time zone to be used in output
#'
#' @return a list with fields sunrise time, sunset time, day length, night
#'   length. The times are returned in the same TZ as used for the date.
#'
#' @keywords manip misc
#' @export
#' @examples
#' library(lubridate)
#' day_night()
#' day_night(ymd("2014-05-30"), lat = 30, lon = 0)
#' day_night(ymd("2014-05-30"), lat = 30, lon = 0, twilight = "civil")
#' day_night(ymd("2014-05-30"), lat = 30, lon = 0, twilight = -6)

day_night <- function(t = today(), lon = 0, lat = 0, twilight = "none", tz=Sys.timezone()) {
  if (!is.numeric(twilight)) {
    if (twilight=="none") {
      twilight_angle <- 0
    } else if (twilight=="civil") {
      twilight_angle <- -6
    } else if (twilight=="nautical") {
      twilight_angle <- -12
    } else if (twilight=="astronomical") {
      twilight_angle <- -18
    } else {
      twilight_angle <- NA
    }
  } else {
    twilight_angle <- ifelse(twilight < 90 & twilight > -33, twilight, NA)
  }
  if (is.na(twilight_angle)) {
    warning("Unrecognized argument for 'twilight' :", twilight)
    return(NA)
  }
  if (!is.POSIXct(t)) {
    if (is.instant(t)) {
      t <- as.POSIXct(t, tz="UTC")
    } else {
      warning("t is not a valid time or date")
    }
  }
  t <- as.POSIXct(t, tz=tz)
  hour(t) <- 0
  minute(t) <- 0
  second(t) <- 0
  t_num <- as.numeric(t, tz=tz)

  altitude <- function(x){
    t_temp <- as.POSIXct(x, origin=origin, tz="UTC")
    return(sun_angles(t_temp,
                      lon=lon,
                      lat=lat)$elevation -
             twilight_angle)
  }

  noon <- try(
    optimize(f=altitude, interval=c(t_num + 7200, t_num + 86400 - 7200), maximum=TRUE)$maximum
  )
  if (inherits(noon, "try-error")) {
    return(NA)
  }
  noon_time <- as.POSIXct(noon, tz=tz, origin=origin)

  rise <- try(
    uniroot(f=altitude, lower = noon - 86400/2, upper=noon)$root,
    silent=TRUE)
  if (inherits(rise, "try-error")) {
    rise <- NA # never
  }
  rise_time <- as.POSIXct(rise, tz=tz, origin=origin)

  set <- try(
    uniroot(f=altitude, lower=noon, upper=noon + 86400/2)$root,
    silent=TRUE)
  if (inherits(set, "try-error")) {
    set <- NA # never
  }
  set_time <- as.POSIXct(set, tz=tz, origin=origin)

  if (is.na(rise) || is.na(set)) {
    daylength <- ifelse(altitude(noon) > 0, 24, 0)
  } else {
    daylength <- set_time - rise_time
  }
return(list(day         = t,
            sunrise     = rise_time,
            noon        = noon_time,
            sunset      = set_time,
            daylength   = daylength,
            nightlength = 24 - daylength
))
}
