# This are just benchmarks
# We need to check also if values are correct
library(tibble)
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggpmisc)
library(photobiology)
library(solartime)
library(suncalc)
library(microbenchmark)
library(fishmethods)
library(svglite)

## vectorised time points

num.years <- 1
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
  num.evals <- as.integer(max(7, min(21, length(times.per.hour) / length(t))))
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
  expand_limits(y = 0) +
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
t <- as.POSIXct(seq(from = today(tzone = "EET"), to = today(tzone = "EET") + days(1),
                    length.out =  24))

benchmark.results <- list()
idx <- 0L
for (l in latitudes) {
  idx <- idx + 1L
  num.evals <- as.integer(max(11, min(99, 5 * length(latitudes[[1]]) / length(l))))
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
result.vec.latitudes.tb %>%
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
  expand_limits(y = 0) +
  theme_bw() + theme(legend.position = "top") -> per.timepoint.latitudes.ggp

svglite("./test/per-timepoint-latitudes-ggp.svg", width = 6, height = 4)
print(per.timepoint.latitudes.ggp)
dev.off()

## profile

t <- as.POSIXct(seq(from = today(tzone = "EET"), to = today(tzone = "EET") + days(1),
                    length.out = 3))
localities <- data.frame(lat = latitudes$thousands, lon = 0)
profvis::profvis(sun_angles(t, geocode = localities))

## test values computed for extreme dates

date_times <-  as.POSIXct((ymd_hm("0001-01-01 00:00") +
                             years(c(0, 100, 1000, 1800, 1900, 2000, 2020, 2100, 3000) - 1)) +
                            days(rep((0:11) * 30, 9)) +
                            hours(rep(c(0,6,12,18), 9 * 12)))
localities <- data.frame(lat = latitudes$tens, lon = 0)[-1, ] # computeSunPosition fails at the equator

length(date_times) * nrow(localities)

## this need to be expanded to multiple times and locations
time <- date_times[7] + months(2)
astrocalc4r(day = day(time), month = month(time), year = year(time), hour = hour(time) + 0.06,
             timezone = 0, lat = localities[50, "lat"], localities[50, "lon"],
            withinput = TRUE, seaorland = "continental")

day_night(time + minutes(3) + seconds(40), geocode = localities[50, ])

sun_angles(time + minutes(3) + seconds(40), geocode = localities[50, ])

getSunlightPosition(time + minutes(3) + seconds(40),
                    lat = localities[50, "lat"], lon = localities[50, "lon"]) %>%
  mutate(altitude = altitude * 180 / pi,
         azimuth = (pi + azimuth) * 180 / pi %% 360)

## photobiology

sa_results.tb <- cbind(sun_angles(date_times, geocode = localities),
                           dplyr::select(day_night(date_times, geocode = localities),
                                  -day, -tz, -longitude, -latitude, -address))
names(sa_results.tb) <- paste("sa.", names(sa_results.tb), sep = "")

nrow(sa_results.tb)

###

csp_results.lst <- list(length(localities))
sum_nrows <- 0
for (i in 1:nrow(localities)) {
  with(localities[i, ],
       computeSunPosition(date_times, lat, lon)) %>%
    as.data.frame() %>%
    mutate(declination = declination * 180 / pi,
           elevation = elevation * 180 / pi,
           azimuth = azimuth * 180 / pi) -> tmp
  sum_nrows <- sum_nrows + nrow(tmp)
  csp_results.lst[[i]] <- tmp
}
print(sum_nrows)

csp_results.tb <- dplyr::bind_rows(csp_results.lst)
names(csp_results.tb) <- paste("csp.", names(csp_results.tb), sep = "")

##

sc_results.lst <- list(length(localities))
sum_nrows <- 0
for (i in 1:nrow(localities)) {
  with(localities[i, ],
       getSunlightPosition(date_times,
                           lat = lat, lon = lon)) %>%
         mutate(altitude = altitude * 180 / pi,
                azimuth = (pi + azimuth) * 180 / pi %% 360) -> tmp
  sum_nrows <- sum_nrows + nrow(tmp)
  sc_results.lst[[i]] <- tmp
}
print(sum_nrows)

sc_results.tb <- dplyr::bind_rows(sc_results.lst)
names(sc_results.tb) <- paste("sc.", names(sc_results.tb), sep = "")

##

ac4r_results.lst <- list(length(localities) * length(date_times))
sum_nrows <- 0
for (i in 1:nrow(localities)) {
  for (j in 1:length(date_times)) {
    astrocalc4r(day = day(date_times[j]), month = month(date_times[j]), year = year(date_times[j]), hour = hour(date_times[j]),
                timezone = 0, lat = localities[["lat"]][1], lon = localities[["lon"]][1],
                withinput = TRUE, seaorland = "continental") -> tmp
  sum_nrows <- sum_nrows + nrow(tmp)
  ac4r_results.lst[[(i - 1) * length(date_times) + j]] <- tmp
  }
}
print(sum_nrows)

ac4r_results.tb <- dplyr::bind_rows(ac4r_results.lst)
names(ac4r_results.tb) <- paste("acr.", names(ac4r_results.tb), sep = "")

##

all_results.tb <- bind_cols(sa_results.tb, csp_results.tb, sc_results.tb, ac4r_results.tb) %>%
  mutate(sa_year = year(sa.time),
         sa_hour = hour(sa.time),
         sa_hour_solar = round(sa.solartime),
         sa_hour_solar = ifelse(sa_hour_solar == 24, 0, sa_hour_solar),csp_delta_elevation = csp.elevation - sa.elevation,
         csp_delta_azimuth = csp.azimuth - sa.azimuth,
         # csp_delta_azimuth = ifelse(csp_delta_azimuth < -180, NA, csp_delta_azimuth),
         # csp_delta_azimuth = ifelse(csp_delta_azimuth >  180, NA, csp_delta_azimuth),
         csp_delta_elevation = csp.elevation - sa.elevation,
         csp.solartime = csp.hour,
         csp_delta_solar_time = ifelse(csp.solartime < 0, 24 - csp.solartime, csp.solartime) - sa.solartime,

         sc_delta_azimuth = sc.azimuth - acr.azimuth,
         sc_delta_elevation = sc.altitude - sa.elevation,

         acr_delta_azimuth = acr.azimuth - sa.azimuth,
         acr.elevation = -(acr.zenith - 90),
         acr_delta_elevation = acr.elevation - sa.elevation,
         acr_delta_declination = acr.declin - sa.declination,
         acr_delta_daylength = acr.daylight - sa.daylength,
         acr_delta_sunset = acr.sunset - sa.sunset,
         acr_delta_sunrise = acr.sunrise - sa.sunrise,
         acr_delta_noon = acr.noon - sa.noon
  )

### Azimuth

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, acr_delta_azimuth, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar azimuth (degrees)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> acr4r_azimuth_error.ggp

svglite("./test/acr4r-azimuth-error-ggp.svg", width = 6, height = 10)
print(acr4r_azimuth_error.ggp)
dev.off()

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, csp_delta_azimuth, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar azimuth (degrees)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: computeSunPosition() vs. sun_angles()") -> csp_azimuth_error.ggp

svglite("./test/csp-azimuth-error-ggp.svg", width = 6, height = 10)
print(csp_azimuth_error.ggp)
dev.off()


all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, sc_delta_azimuth, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar azimuth (degrees)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: getSunlightPosition() vs. sun_angles()") -> sc_azimuth_error.ggp

svglite("./test/sc-azimuth-error-ggp.svg", width = 6, height = 10)
print(sc_azimuth_error.ggp)
dev.off()

### Elevation

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, acr_delta_elevation, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar elevation (degrees)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> acr4r_elevation_error.ggp

svglite("./test/acr4r-elevation-error-ggp.svg", width = 6, height = 10)
print(acr4r_elevation_error.ggp)
dev.off()

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, csp_delta_elevation, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar elevation (degrees)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: computeSunPosition() vs. sun_angles()") -> csp_elevation_error.ggp

svglite("./test/csp-elevation-error-ggp.svg", width = 6, height = 10)
print(csp_elevation_error.ggp)
dev.off()


all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, sc_delta_elevation, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar elevation (degrees)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: getSunlightPosition() vs. sun_angles()") -> sc_elevation_error.ggp

svglite("./test/sc-elevation-error-ggp.svg", width = 6, height = 10)
print(sc_elevation_error.ggp)
dev.off()


#### Solar time

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, csp_delta_solar_time * 60, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in solar time (min)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: computeSunPosition() vs. sun_angles()") -> acr_solartime_error.ggp

svglite("./test/acr-solartime-error-ggp.svg", width = 6, height = 10)
print(acr_solartime_error.ggp)
dev.off()

### Sunrise

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, acr_delta_sunrise * 60, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in sunrise time (min)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> acr4r_sunrise_error.ggp

svglite("./test/acr4r-sunrise-error-ggp.svg", width = 6, height = 10)
print(acr4r_sunrise_error.ggp)
dev.off()

### Noon

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, acr_delta_noon * 60, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in noon time (min)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> acr4r_noon_error.ggp

svglite("./test/acr4r-noon-error-ggp.svg", width = 6, height = 10)
print(acr4r_noon_error.ggp)
dev.off()

### Sunset

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, acr_delta_sunset * 60, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in sunset time (min)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> acr4r_sunset_error.ggp

svglite("./test/acr4r-sunset-error-ggp.svg", width = 6, height = 10)
print(acr4r_sunset_error.ggp)
dev.off()

### Daylength

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, acr_delta_daylength * 60, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(alpha = 1/3) +
  stat_quadrant_counts(size = 3) +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Difference in daylength time (min)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> acr4r_daylength_error.ggp

svglite("./test/acr4r-daylength-error-ggp.svg", width = 6, height = 10)
print(acr4r_daylength_error.ggp)
dev.off()

### Actual values

all_results.tb %>%
  filter(sa_year != 3000 & sa_year >= 1000) %>%
  ggplot(aes(sa.latitude, colour = factor(sa_hour))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(mapping = aes(y = sa.sunset), alpha = 1/3) +
  geom_point(mapping = aes(y = acr.sunset), alpha = 1/3, shape = "triangle") +
  expand_limits(y = c(-2, 2)) +
  facet_wrap(~sa_year, labeller = label_both, ncol = 2) +
  labs(y = "Sunset time (h)",
       x = "Latitude (degrees)",
       colour = "Time of day\nUTC (h)") +
  theme_bw() + theme(legend.position = "top") +
  ggtitle("Azimuth: astrocalc4r() vs. sun_angles()") -> sunset.ggp

svglite("./test/sunset-ggp.svg", width = 6, height = 10)
print(sunset.ggp)
dev.off()

## other functions

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
