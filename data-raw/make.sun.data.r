library(photobiology)
library(dplyr)
setwd("data-raw")
sun.midday.data <- read.table("sun_20100622_midday.txt", col.names=c("time_min","w.length","s.e.irrad"))
sun.midday.data$time_min <- NULL
sun.midday.data$w.length <- sun.midday.data$w.length / 10.0
sun.midday.data$s.e.irrad <- sun.midday.data$s.e.irrad / 1e3
sun.midday.data$s.q.irrad <- with(sun.midday.data, as_quantum_mol(w.length, s.e.irrad))

sun.data <- as_data_frame(sun.midday.data)
sun.spct <- as.source_spct(sun.spct)

setwd("../data")

save(sun.spct, file="sun.spct.rda")
save(sun.data, file="sun.data.rda")
rm(sun.data, sun.spct, sun.midday.data)

setwd("..")

setwd("data-raw")

sun.daily.data <- read.table("sun_20120601_cum.hel.txt", col.names=c("w.length","s.e.irrad"))
sun.daily.data$w.length <- sun.daily.data$w.length / 10.0
# sun.daily.data$s.e.irrad <- sun.daily.data$s.e.irrad / 1e3
sun.daily.data$s.q.irrad <- with(sun.daily.data, as_quantum_mol(w.length, s.e.irrad))
sun.daily.data <- as_data_frame(sun.daily.data)
sun.daily.spct <- as.source_spct(sun.daily.data, time.unit = "day")

setwd("../data")

save(sun.daily.spct, file="sun.daily.spct.rda")
save(sun.daily.data, file="sun.daily.data.rda")
rm(sun.daily.spct, sun.daily.data)

setwd("..")
