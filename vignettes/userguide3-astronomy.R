## ---- include=FALSE, echo=FALSE------------------------------------------
knitr::opts_chunk$set(fig.width=8, fig.height=4)

## ---- printing-spectra, eval=TRUE, include=FALSE-------------------------
# library(tibble)
options(tibble.print_max = 6, tibble.print_min = 4)

## ---- pkg-load, eval=TRUE------------------------------------------------
library(photobiology)
library(lubridate)

## ------------------------------------------------------------------------
sun_angles(now(), geocode = data.frame(lat = 34, lon = 0))
sun_angles(ymd_hms("2014-01-01 0:0:0", tz = "UTC") + hours(1:3))

## ------------------------------------------------------------------------
sun_angles(getWhenMeasured(sun.spct), geocode = getWhereMeasured(sun.spct))

## ------------------------------------------------------------------------
sun_elevation(ymd_hms("2014-01-01 0:0:0", tz = "UTC") + hours(1:3))

## ------------------------------------------------------------------------
sun_zenith_angle(ymd_hms("2014-01-01 0:0:0", tz = "UTC") + hours(1:3))

## ------------------------------------------------------------------------
sun_azimuth(ymd_hms("2014-01-01 0:0:0", tz = "UTC") + hours(1:3))

## ------------------------------------------------------------------------
dates <- seq(from = ymd("2015-03-01"), to = ymd("2015-07-1"), length.out = 3)

## ------------------------------------------------------------------------
noon_time(dates, tz = "UTC", data.frame(lat = 34, lon = 0))

## ------------------------------------------------------------------------
noon_time(dates, tz = "CET", data.frame(lat = 34, lon = 0))

## ------------------------------------------------------------------------
day_night(dates, geocode = data.frame(lat = 60, lon = 0))

## ------------------------------------------------------------------------
sunrise_time(geocode = data.frame(lat = 60, lon = 0))

## ------------------------------------------------------------------------
sunrise_time(today("UTC"), tz = "UTC", geocode = data.frame(lat = 60, lon = 0))
sunrise_time(today("EET"), tz = "EET", geocode = data.frame(lat = 60, lon = 25))

## ------------------------------------------------------------------------
sunrise_time(dates, geocode = data.frame(lat = 60, lon = 0))
sunrise_time(dates, geocode = data.frame(lat = -60, lon = 0))

## ------------------------------------------------------------------------
sunrise_time(today("EET"), tz = "EET", 
             geocode = data.frame(lat = 60, lon = 25),
             twilight = "civil")
sunrise_time(today("EET"), tz = "EET", 
             geocode = data.frame(lat = 60, lon = 25),
             twilight = -10)
sunrise_time(today("EET"), tz = "EET", 
             geocode = data.frame(lat = 60, lon = 25),
             twilight = +12)

## ------------------------------------------------------------------------
sunrise_time(today("EET"), 
             tz = "EET", 
             geocode = data.frame(lat = 60, lon = 25),
             unit.out = "hours")

## ------------------------------------------------------------------------
day_length(dates, geocode = data.frame(lat = 60, lon = 25))
night_length(dates, geocode = data.frame(lat = 60, lon = 25))

## ------------------------------------------------------------------------
day_night(dates, 
          geocode = data.frame(lat = 60, lon = 25))
day_night(dates, 
          geocode = data.frame(lat = 60, lon = 25), 
          unit.out = "datetime")

## ------------------------------------------------------------------------
Paris.geo <- data.frame(lon = 2.352222, lat = 48.85661, address = "Paris")
Paris.time <- ymd_hms("2016-09-30 06:00:00", tz = "UTC")
solar_time(Paris.time, geocode = Paris.geo)
solar_time(Paris.time, geocode = Paris.geo, unit.out = "datetime")

## ------------------------------------------------------------------------
my.solar.t <- solar_time(Paris.time, geocode = Paris.geo)
is.solar_time(my.solar.t)
is.numeric(my.solar.t)

## ------------------------------------------------------------------------
my.solar.d <- solar_time(Paris.time, geocode = Paris.geo, unit.out = "datetime")
is.solar_date(my.solar.d)
is.timepoint(my.solar.d)

## ------------------------------------------------------------------------
times <- now() + days(0:1)
times
as_tod(times)
as_tod(times, unit.out = "minutes")

