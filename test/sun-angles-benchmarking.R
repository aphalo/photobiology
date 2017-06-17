## test and benchmark sun_angles()

library(lubridate)
library(photobiology)
library(microbenchmark)
library(profvis)

str(sun_angles())
str(sun_angles(now() + minutes(1:5)))

microbenchmark(sun_angles(now() + minutes(1:1000)))
profvis(sun_angles(now() + minutes(1:1000)))

