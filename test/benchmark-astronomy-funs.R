# This are just benchmarks
# We need to check also if values are correct
library(tibble)
library(dplyr)
library(lubridate)
library(ggplot2)
library(photobiology)
library(solartime)
library(microbenchmark)
library(svglite)

## vectorised time points

num.years <- 10
times.per.min <- as.POSIXct(seq(from = today(tzone = "EET"), to = today(tzone = "EET") + years(num.years),
                 length.out =  365 * 24 * 60 * (num.years)))
times.per.hour <- as.POSIXct(seq(from = today(tzone = "EET"), to = today(tzone = "EET") + years(num.years),
                             length.out =  365 * 24 * (num.years)))
times.per.day <- as.POSIXct(seq(from = today(tzone = "EET"), to = today(tzone = "EET") + years(num.years),
                             length.out =  365 * (num.years)))
times.per.month <- as.POSIXct(seq(from = today(tzone = "EET"), to = today(tzone = "EET") + years(num.years),
                            length.out =  12 * (num.years)))
times.per.year <- as.POSIXct(seq(from = today(tzone = "EET"), to = today(tzone = "EET") + years(num.years),
                             length.out = num.years))
length(times.per.year)
times.per.year

times <- list(per.min = times.per.min, per.hour = times.per.hour, per.day = times.per.day, per.month = times.per.month, per.year = times.per.year)

benchmark.results <- list()
idx <- 0L
for (t in times) {
  idx <- idx + 1L
  num.evals <- as.integer(max(7, min(99, length(times.per.hour) / length(t))))
  result <-
    microbenchmark(sun_angles(t, geocode = data.frame(lat = 60, lon = 0)),
                   computeSunPosition(timestamp = t, 60, 0),
                   times = num.evals, unit = "s",
                   setup = gc())
  result$case <- names(times)[idx]
  result$case.length <- length(t)
  benchmark.results <- c(benchmark.results, list(as.data.frame(result)))
  print(result)
    cat(".")
}
result.vec.times.tb <- bind_rows(benchmark.results)
bind_rows(benchmark.results) %>%
  mutate(fun = ifelse(grepl("sun_angles", expr), "photobiology::sun_angles", "solartime::computeSunPosition")) %>%
  group_by(case, fun) %>%
  summarise(case.length = median(case.length),
            median.time = median(time) * 1e-9, # seconds
            median.time.per.time.point = median.time / case.length) -> summary.vec.times.tb
summary.vec.times.tb

summary.vec.times.tb %>%
  ggplot(aes(case.length, median.time * 1000, colour = fun)) +
  geom_point() +
  geom_line() +
  scale_x_log10(name = "Length of vector of times points") +
  scale_y_log10(name = "Total execution time (ms)") +
  theme_bw() + theme(legend.position = "top") -> total.times.ggp

svglite("./test/total-times-ggp.svg", width = 6, height = 4)
print(total.times.ggp)
dev.off()

summary.vec.times.tb %>%
  ggplot(aes(case.length, median.time.per.time.point * 1e6, colour = fun)) +
  geom_point() +
  geom_line() +
  scale_x_log10(name = "Length of vector of time points") +
  scale_y_continuous(name = expression("Execution time per time point"~~(mu*s))) +
  theme_bw() + theme(legend.position = "top") -> per.timepoint.ggp

svglite("./test/per-timepoint-ggp.svg", width = 6, height = 4)
print(per.timepoint.ggp)
dev.off()

## Vectorised latitudes

latitudes <- list(thousands = seq(0, 87, length.out = 2^12),
                  hundreds = seq(0, 87, length.out = 2^9),
                  tens = seq(0, 87, length.out = 2^6),
                  ones = seq(0, 87, length.out = 2^3))

# every hour for a week
t <- as.POSIXct(seq(from = today(tzone = "EET"), to = today(tzone = "EET") + weeks(1),
                    length.out =  7 * 24))

benchmark.results <- list()
idx <- 0L
for (l in latitudes) {
  idx <- idx + 1L
  num.evals <- as.integer(max(7, min(99, length(latitudes$thousands) / length(t))))
  localities <- data.frame(lat = l, lon = 0)
  result <-
    microbenchmark(sun_angles(t, geocode = localities),
                   for (i in 1:nrow(localities)) {
                     with(localities[i, ],
                          computeSunPosition(timestamp = t, lat, lon)
                     )
                   },
                   times = num.evals, unit = "s",
                   setup = gc())
  result$case <- names(latitudes)[idx]
  result$case.length <- length(l) * length(t)
  result$num.localities <- length(l)
  benchmark.results <- c(benchmark.results, list(as.data.frame(result)))
  print(result)
  cat(".")
}
result.vec.latitudes.tb <- bind_rows(benchmark.results)
bind_rows(benchmark.results) %>%
  mutate(fun = ifelse(grepl("sun_angles", expr), "photobiology::sun_angles", "solartime::computeSunPosition")) %>%
  group_by(case, fun) %>%
  summarise(case.length = median(case.length),
            num.localities = median(num.localities),
            median.time = median(time) * 1e-9, # seconds
            median.time.per.time.point = median.time / case.length) -> summary.vec.latitudes.tb
summary.vec.latitudes.tb

summary.vec.latitudes.tb %>%
  ggplot(aes(num.localities, median.time * 1000, colour = fun)) +
  geom_point() +
  geom_line() +
  scale_x_log10(name = "Length of vector of latitudes") +
  scale_y_log10(name = "Total execution time (ms)") +
  theme_bw() + theme(legend.position = "top") -> total.times.latitudes.ggp

svglite("./test/total-times-latitudes-ggp.svg", width = 6, height = 4)
print(total.times.latitudes.ggp)
dev.off()

summary.vec.latitudes.tb %>%
  ggplot(aes(num.localities, median.time.per.time.point * 1e6, colour = fun)) +
  geom_point() +
  geom_line() +
  scale_x_log10(name = "Length of vector of latitudes") +
  scale_y_continuous(name = expression("Execution time per time point"~~(mu*s))) +
  theme_bw() + theme(legend.position = "top") -> per.timepoint.latitudes.ggp

svglite("./test/per-timepoint-latitudes-ggp.svg", width = 6, height = 4)
print(per.timepoint.latitudes.ggp)
dev.off()

## test values computed for extreme dates

date_times <-  as.POSIXct(ymd_hm("2019-05-01 15:30") + years(c(-4000L, -1000L, -100, 0, 100, 1000)))

sun_angles(date_times, geocode = data.frame(lat = 5, lon = 0))
computeSunPosition(date_times, 5, 0) %>%
  as.data.frame() %>%
  mutate(declination = declination * 180 / pi,
         elevation = elevation * 180 / pi,
         azimuth = azimuth * 180 / pi)

microbenchmark(sun_angles_photobiology.df <- sun_angles(times),
               times = 30, unit = "s")
sun_angles_photobiology.df <- sun_angles(times)

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
