## test and benchmark sun_angles()

library(tidyverse)
library(lubridate)
library(photobiology)
library(microbenchmark)
library(profvis)

str(sun_angles())
str(sun_angles(now() + minutes(1:5)))

microbenchmark(sun_angles(now("UTC") + minutes(1:1000)) -> z, unit = "ms")
z

times <- data.frame(size = c(1, 2, 10, 100, 500, 1000, 5000, 10000, 50000),
           ms = c(3.77, 3.80, 3.84, 4.42, 6.54, 9.63, 29.81, 57.2, 271))

ggplot(times, aes(size, ms)) +
  geom_line()

profvis(sun_angles(now() + minutes(1:10000)))

my.geocodes <- data.frame(lat = c(0, 15, 30, 45, 60, 89.99),
                          lon = 0,
                          address = paste("latitude",  c(0, 15, 30, 45, 60, 89.99), sep = "."))

my.lon.geocodes <- rbind(my.geocodes, my.geocodes)
my.lon.geocodes[["lon"]] <- rep(c(0, 30), rep(6, 2))

microbenchmark(sun_angles(now("UTC") + minutes(1:10), geocode = my.geocodes) -> z, unit = "ms")
z

microbenchmark(sun_angles(now("UTC") + minutes(1:10), geocode = my.lon.geocodes) -> z, unit = "ms")
z

microbenchmark(day_night())
microbenchmark(day_night(geocode = my.geocodes))
microbenchmark(day_night(now(tz = "UTC") + weeks(1:5)))
microbenchmark(day_night(now() + days(1:365), geocode = my.geocodes) -> z, unit = "ms")
z

profvis(day_night(now() + days(1:365), geocode = my.geocodes) -> z)
