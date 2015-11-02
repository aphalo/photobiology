library(photobiology)
library(dplyr)
library(lubridate)
library(ggmap)

setwd("data-raw")
sun.midday.data <- read.table("sun_20100622_midday.txt", col.names = c("time_min","w.length","s.e.irrad"))
sun.midday.data$time_min <- NULL
sun.midday.data$w.length <- sun.midday.data$w.length / 10.0
sun.midday.data$s.e.irrad <- sun.midday.data$s.e.irrad / 1e3
sun.midday.data$s.q.irrad <- with(sun.midday.data, as_quantum_mol(w.length, s.e.irrad))

sun.data <- as_data_frame(sun.midday.data)
sun.spct <- as.source_spct(sun.data)
sun.spct <- trim_spct(sun.spct, low.limit = 280, fill = 0, use.hinges = FALSE)
setWhenMeasured(sun.spct, ymd_hms("2010-06-22 9:51:00"), tz = "UTC")
setWhereMeasured(sun.spct, geocode("Kumpula, Helsinki, FI", "latlona"))
comment(sun.spct) <- "Simulated solar spectrum based on real weather conditions.\nData author: Dr. Anders Lindfors\nFinnish Meteorological Institute.\nSee help file for references."

setwd("../data")

save(sun.spct, file = "sun.spct.rda")
save(sun.data, file = "sun.data.rda")
rm(sun.data, sun.spct, sun.midday.data)

setwd("..")

setwd("data-raw")

sun.daily.data <- read.table("sun_20120601_cum.hel.txt", col.names = c("w.length","s.e.irrad"))
sun.daily.data$w.length <- sun.daily.data$w.length / 10.0
sun.daily.data$s.q.irrad <- with(sun.daily.data, as_quantum_mol(w.length, s.e.irrad))
sun.daily.data <- as_data_frame(sun.daily.data)
sun.daily.spct <- as.source_spct(sun.daily.data, time.unit = "day")
sun.daily.spct <- trim_spct(sun.daily.spct, low.limit = 280, fill = 0, use.hinges = FALSE)
setWhenMeasured(sun.daily.spct, ymd("2012-06-01"), tz = "UTC")
setWhereMeasured(sun.daily.spct, geocode("Kumpula, Helsinki, FI", "latlona"))
comment(sun.daily.spct) <- "Total daily spectral exposure estimated from hourly simulations of the solar spectrum\nbased real weather conditions.\nData author: Dr. Anders Lindfors\nFinnish Meteorological Institute.\nSee help file for references."

setwd("../data")

save(sun.daily.spct, file = "sun.daily.spct.rda")
save(sun.daily.data, file = "sun.daily.data.rda")
rm(sun.daily.spct, sun.daily.data)

setwd("..")
