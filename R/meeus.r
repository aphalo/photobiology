#' @rdname julian_day
#'
calendar_change <- lubridate::ymd_hms("1582-10-14 22:20:11 UTC") # Julian -> Gregorian

# All functions defined in this file are "internal" and not exported
# They are organized as very small functions to allow reuse of the results of
# partial calculations. All constants are contained in the code itself.
# They are based in "NOAA Sunrise/Sunset and Solar Position Calculators"
# available at http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
#' Solar astronomy using Meeus' algorithm
#'
#' The exact julian day computation is adapted from ode::julianDay() and tested
#' againts test cases in Redas and Andreas (2008, table A4.1) for validity up
#' to year 4712 BC.
#'
#' Low level functions based on NOAA's Excel worksheet
#'
#' @param time dateTime
#' @param x numeric Julian century
#' @param anom numeric Solar anomaly in degrees
#' @param eccent numeric Eccentricity of Earth orbit
#' @param eclip numeric Ecliptic
#' @param app.lon,obliq.corr,mean.lon,nag,decline numeric Angles in degrees
#' @param eq.of.time,ha.sunrise,noon numeric
#' @param zenith.angle,elevation.angle,hour.angle numeric Angles in degrees
#' @param lat,lon numeric Geographic coordinates in degrees
#'
#' @keywords internal
#'
julian_day <- function(time) {
  time <- lubridate::with_tz(time, tzone = "UTC")
  y <- as.double(lubridate::year(time))
  m <- as.double(lubridate::month(time))
  d <- as.double(lubridate::day(time))
  selector <- m <= 2
  m[selector] <- m[selector] + 12
  y[selector] <- y[selector] - 1
  d <- d + as_tod(time, unit.out = "hours", tz = "UTC") / 24
  a <- floor(y / 100)
  b <- 2 - a + floor(a / 4)
  jd <- floor(365.25 * y) + floor(30.6001 * (m + 1)) + d + 1720994.5
  ifelse(time > calendar_change, jd + b, jd)
}

#' @rdname julian_day
#'
julian_day_fast <- function(time) {
  2440587.79166667 + as.numeric(julian(time))
}

#' @rdname julian_day
#'
julian_century <- function(time) {
  (julian_day(time) - 2451545) / 36525
}

#' @rdname julian_day
#'
geom_mean_lon_sun <- function(x) {
  (280.46646 + x * (36000.76983 + x * 0.0003032)) %% 360
}

#' @rdname julian_day
#'
geom_mean_anom_sun <- function(x) {
  357.52911 + x * (35999.05029 - 0.0001537 * x)
}

#' @rdname julian_day
#'
eccent_earth_orbit <- function(x) {
  0.016708634 - x * (0.000042037 + 0.0000001267 * x)
}

#' @rdname julian_day
#'
sun_eq_of_ctr <- function(x, anom) {
  anom.rad <- anom / 180 * pi
  sin(anom.rad) * (1.914602 - x * (0.004817 + 0.000014 * x )) +
    sin(2 * anom.rad) * (0.019993 - 0.000101 * x) +
    sin(3 * anom.rad) * 0.000289
}

#' @rdname julian_day
#'
sun_rad_vector <- function(eccent, anom) {
  anom.rad <- anom / 180 * pi
  (1.000001018*  (1 - eccent^2))/(1 + eccent * cos(anom.rad))
}

#' @rdname julian_day
#'
sun_app_lon <- function(x, lon) {
  lon - 0.00569 - 0.00478 * sin((125.04 - 1934.136 * x) / 180 * pi)
}

#' @rdname julian_day
#'
mean_obliq_eclip <- function(x) {
  23 + (26 + ((21.448 -
                 x * (46.815 + x * (0.00059 - x * 0.001813)))) / 60) / 60
}

#' @rdname julian_day
#'
obliq_corr <- function(x, eclip) {
  eclip + 0.00256 * cos((125.04 - 1934.136 * x) / 180 * pi)
}

#' @rdname julian_day
#'
sun_rt_ascen <- function(app.lon, obliq.corr) {
  app.lon.rad <- app.lon / 180 * pi
  obliq.corr.rad <- obliq.corr / 180 * pi
  atan2(cos(obliq.corr.rad) * sin(app.lon.rad),
        cos(app.lon.rad )) / pi * 180
}

#' @rdname julian_day
#'
sun_decline <- function(app.lon, obliq.corr) {
  app.lon.rad <- app.lon / 180 * pi
  obliq.corr.rad <- obliq.corr / 180 * pi
  asin(sin(obliq.corr.rad) * sin(app.lon.rad)) /
    pi * 180
}

#' @rdname julian_day
#'
var_y <- function(obliq.corr) {
  x <- obliq.corr / 360 * pi
  tan(x)^2
}

#' @rdname julian_day
#'
eq_of_time <- function(mean.lon,
                       eccent.earth,
                       anom.mean,
                       var.y) {
  mean.lon.rad <- mean.lon / 180 * pi
  anom.mean.rad <- anom.mean / 180 * pi
  sin.anom.mean.rad <- sin(anom.mean.rad) # avoid computing it twice
  4 * ((var.y * sin(2 * mean.lon.rad) -
         2 * eccent.earth * sin.anom.mean.rad +
         4 * eccent.earth * var.y * sin.anom.mean.rad * cos(2 * mean.lon.rad) -
         0.5 * var.y^2 * sin(4 * mean.lon.rad) -
         1.25 * eccent.earth^2 * sin(2 * anom.mean.rad)) / pi * 180)
}

#' @rdname julian_day
#'
ha_sunrise <- function(lat, decline, nag = 0) {
  lat.rad <- lat / 180 * pi
  declin.rad <- decline / 180 * pi
  suppressWarnings( # NaNs can be produced for polar regions
    z <-
      acos(cos((90 + nag) / 180 * pi) /            # 90.833 in NOAA's code
             (cos(lat.rad) * cos(declin.rad)) -
             tan(lat.rad) * tan(declin.rad))  / pi * 180)
  z
}

#' @rdname julian_day
#'
solar_noon <- function(lon, eq.of.time) {
  (720 - 4 * lon - eq.of.time) / 1440
}

#' @rdname julian_day
#'
sunrise <- function(noon, ha.sunrise) {
  (noon * 1440 -
     ha.sunrise * 4) / 1440
}

#' @rdname julian_day
#'
sunset <- function(noon, ha.sunrise) {
  (noon * 1440 +
     ha.sunrise * 4) / 1440
}

#' @rdname julian_day
#'
sunlight_duration <- function(ha.sunrise, unit.out = "hours") {
  8 * ha.sunrise
}

#' @rdname julian_day
#'
#' @return datetime
solar_datetime <- function(time, lat, lon, eq.of.time) {
  time + lubridate::seconds((eq.of.time + 4 * lon) * 60)
}

#' @rdname julian_day
#'
#' @return numeric
solar_tod <- function(time, lat, lon, eq.of.time) {
  stopifnot(lubridate::is.instant(time))
  tod <- as_tod(time, unit.out = "minutes", tz = "UTC")
  (tod + eq.of.time + 4 * lon) %% 1440
}

#' @rdname julian_day
#'
hour_angle <- function(solar.time) {
  ifelse(solar.time < 0, solar.time / 4 + 180, solar.time / 4 - 180)
}

#' @rdname julian_day
#'
zenith_angle <- function(lat, hour.angle, decline) {
  lat.rad <- lat / 180 * pi
  hour.angle.rad <- hour.angle / 180 * pi
  declin.rad <- decline / 180 * pi
  (acos(sin(lat.rad) * sin(declin.rad) +
         cos(lat.rad) * cos(declin.rad) * cos(hour.angle.rad))) / pi * 180
}

#' @rdname julian_day
#'
elevation_angle <- function(lat, hour.angle, decline) {
  90 - zenith_angle(lat, hour.angle, decline)
}

#' @rdname julian_day
#'
atm_refraction_approx <- function(elevation.angle) {
# The effects of variation in atmospheric pressure and temperature are ignored
# in this implementation.
#
  elev.rad <- elevation.angle / 180 * pi
  ifelse(elevation.angle > 85,
         0,
         ifelse(elevation.angle > 5,
                58.1 / tan(elev.rad) - 0.07 / tan(elev.rad)^3 + 0.000086 / tan(elev.rad)^5,
                ifelse(elevation.angle > -0.575,
                       1735 + elevation.angle *
                         (-518.2 + elevation.angle * (103.4 + elevation.angle * (-12.79 + elevation.angle * 0.711))),
                       -20.772 / tan(elev.rad)))) / 3600
}

#' @rdname julian_day
#'
azimuth_angle <- function(lat, hour.angle, zenith.angle, decline) {
  lat.rad <- lat / 180 * pi
  zenith.angle.rad <- zenith.angle / 180 * pi
  declin.rad <- decline / 180 * pi
#  hour.angle.rad <- hour.angle / 180 * pi
  ifelse(hour.angle > 0,
         (acos(((sin(lat.rad) * cos(zenith.angle.rad)) -
                  sin(declin.rad)) / (cos(lat.rad) * sin(zenith.angle.rad))) /
            pi * 180 + 180) %% 360,
         (540 - (acos(((sin(lat.rad) * cos(zenith.angle.rad)) -
                        sin(declin.rad)) / (cos(lat.rad) * sin(zenith.angle.rad)))) /
                  pi * 180) %% 360)
}

# IF(AC2>0,MOD(DEGREES(ACOS(((SIN(RADIANS($B$3))*COS(RADIANS(AD2)))-
#                 SIN(RADIANS(T2)))/(COS(RADIANS($B$3))*SIN(RADIANS(AD2)))))+180,360),
#    MOD(540-DEGREES(ACOS(((SIN(RADIANS($B$3))*COS(RADIANS(AD2)))-
#                            SIN(RADIANS(T2)))/(COS(RADIANS($B$3))*SIN(RADIANS(AD2))))), 360))
