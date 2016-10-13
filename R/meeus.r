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
#' @param eq.of.time,ha_sunrise,solar.noon.utc numeric
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
                       vary.y) {
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
  acos(cos((90.833 + ang) / 180 * pi) /
         (cos(lat.rad) * cos(declin.rad))-
         tan(lat.rad) * tan(declin.rad))  / pi * 180
}

#' @rdname julian.day
#'
solar_noon_utc <- function(lon, eq.of.time) {
  (720 - 4 * lon - eq.of.time) / 1440
}

#' @rdname julian.day
#'
sunrise_utc <- function(solar.noon.utc, ha.sunrise) {
  ((720 - 4 * solar.noon.utc - eq.of.time) -
     ha.sunrise * 4) / 1440
}

#' @rdname julian.day
#'
sunset_utc <- function(solar.noon.utc, ha.sunrise) {
  ((720 - 4 * solar.noon.utc - eq.of.time) +
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



