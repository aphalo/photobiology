# This are just benchmarks
# We need to check also if values are correct

library(lubridate)
library(ggplot2)
library(microbenchmark)

times <- as.POSIXct(seq(from = today(tzone = "EET"), to = today(tzone = "EET") + years(1),
             length.out =  365 * 24 * 60))
length(times)

library(photobiology)
microbenchmark(sun_angles_photobiology.df <- sun_angles(times),
               times = 30, unit = "s")
sun_angles_photobiology.df <- sun_angles(times)

library(solartime)
microbenchmark(sun_angles_photobiology.df <- sun_angles(times, geocode = data.frame(lat = 60, lon = 0)),
               times = 30, unit = "s")
sun_angles_solartime.df <- computeSunPosition(timestamp = times, 60, 0)

sun_angles.bench <-
  microbenchmark(sun_angles(times, geocode = data.frame(lat = 60, lon = 0)),
                 computeSunPosition(timestamp = times, 60, 0),
                 times = 30, unit = "s")
autoplot(sun_angles.bench)

days <- seq(from = today(tzone = "EET") - years(1000), to = today(tzone = "EET"),
            length.out =  365 * 1000)
length(days)

microbenchmark(sunrise_times_photobiology.df <- day_night(days, geocode = data.frame(lat = 60, lon = 0)),
               times = 30, unit = "s")
sunrise_times_photobiology.df <- day_night(days, geocode = data.frame(lat = 60, lon = 0))

microbenchmark(sunrise_times_solartime.df <- computeSunriseHour(days, 60, 0),
               times = 30, unit = "s")
sunrise_times_solartime.df <- computeSunriseHour(days, 60, 0)

library(suncalc)

microbenchmark(sunrise_times_suncalc.df <- getSunlightTimes(days[1:(365*10)], 60, 0),
               times = 5, unit = "s")
sunrise_times_suncalc.df <- getSunlightTimes(days, 60, 0)

sun_rise.bench <-
  microbenchmark(getSunlightTimes(days[1:(365*50)], 60, 0),
                 computeSunriseHour(days[1:(365*50)], 60, 0),
                 day_night(days[1:(365*50)], geocode = data.frame(lat = 60, lon = 0), unit.out = "day"),
                 day_night(days[1:(365*50)], geocode = data.frame(lat = 60, lon = 0), unit.out = "datetime"),
                 times = 5, unit = "s")

autoplot(sun_rise.bench)

# test values

geocode <- data.frame(lat = 60.17, lon = 24.94, address = "Helsinki, Finland")
day_night(geocode = geocode, tz = "EET", unit.out = "datetime")
day_night(geocode = geocode, tz = "EET", unit.out = "hour")
local_t <- day_night(geocode = geocode, tz = "EET", unit.out = "datetime")$sunrise
local_t
sol_t <- solar_time(local_t, geocode = geocode)
sol_t
computeSunriseHour(today(tzone = "EET"), longDeg = geocode$lon, latDeg = geocode$lat)
