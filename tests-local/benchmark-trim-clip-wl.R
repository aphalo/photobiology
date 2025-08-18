# This are just benchmarks
# We need to check also if values are correct
library(microbenchmark)
library(photobiology)

my.spct <- q2e(sun.spct, action = "replace")

microbenchmark(trim_wl(my.spct, range = c(400, 700), fill = NULL))
microbenchmark(trim_wl(my.spct, range = c(400, 700), fill = NA))
microbenchmark(trim_wl(my.spct, range = c(400, 700), fill = 0))
microbenchmark(trim_wl(my.spct, range = c(400, 700)))
microbenchmark(trim_wl(my.spct, range = c(400, NA)))
microbenchmark(trim_wl(my.spct, range = c(NA, 700)))

microbenchmark(clip_wl(my.spct, range = c(400, 700), fill = NULL))
microbenchmark(clip_wl(my.spct, range = c(400, 700), fill = NA))
microbenchmark(clip_wl(my.spct, range = c(400, 700), fill = 0))
microbenchmark(clip_wl(my.spct, range = c(400, 700)))
microbenchmark(clip_wl(my.spct, range = c(400, NA)))
microbenchmark(clip_wl(my.spct, range = c(NA, 700)))

for (i in 1:1000) trim_wl(my.spct, range = c(400, 700))

