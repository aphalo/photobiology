## test and benchmark sun_angles()

library(tidyverse)
library(tibble)
library(lubridate)
library(photobiology)
library(microbenchmark)
library(profvis)

str(sun_angles())
str(sun_angles(now() + minutes(1:5)))

my.times <- now("UTC") + minutes(1:50000)
microbenchmark(sun_angles(my.times) -> z, unit = "ms")
z

times <- tibble(size = c(1, 2, 10, 100, 500, 1000, 5000, 10000, 50000),
           ms = c(3.77, 3.80, 3.84, 4.42, 6.54, 9.63, 29.81, 57.2, 271))

ggplot(times, aes(size, ms)) +
  geom_line()

my.times <- now("UTC") + minutes(1:5000000)

profvis(sun_angles(my.times))

my.geocodes <- tibble(lat = c(0:89, 89.99),
                          lon = 0,
                          address = paste("latitude",  c(0:89, 89.99), sep = "."))

my.lon.geocodes <- rbind(my.geocodes, my.geocodes)
my.lon.geocodes[["lon"]] <- rep(c(0, 30), rep(91, 2))

my.times <- now("UTC") + minutes(1:10)

microbenchmark(sun_angles(my.times, geocode = my.geocodes) -> z, unit = "ms")
z

profvis(sun_angles(my.times, geocode = my.geocodes))

microbenchmark(sun_angles(my.times, geocode = my.lon.geocodes) -> z, unit = "ms")
z

profvis(sun_angles(my.times, geocode = my.lon.geocodes))

microbenchmark(day_night())
microbenchmark(day_night(geocode = my.geocodes))
microbenchmark(day_night(geocode = my.geocodes, unit.out = "datetime"))

my.times <- now(tzone = "UTC") + weeks(1:5)
microbenchmark(day_night(my.times))

my.times <- now() + days(1:365)
microbenchmark(day_night(my.times, geocode = my.geocodes[1:2, ]) -> z, unit = "ms")
z

my.times <- now() + days(1:365)
profvis(day_night(my.times, geocode = my.geocodes) -> z)
