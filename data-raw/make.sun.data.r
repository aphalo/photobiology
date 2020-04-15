library(photobiology)
library(dplyr)
library(lubridate)

geocode.kumpula <- data.frame(lat = 60.20911,
                              lon = 24.96474,
                              address = "Kumpula, Helsinki, FI",
                              stringsAsFactors = FALSE)

sun.midday.data <- read.table("data-raw/sun_20100622_midday.txt", col.names = c("time_min","w.length","s.e.irrad"))
sun.midday.data$time_min <- NULL
sun.midday.data$w.length <- sun.midday.data$w.length / 10.0
sun.midday.data$s.e.irrad <- sun.midday.data$s.e.irrad / 1e3
sun.midday.data$s.q.irrad <- with(sun.midday.data, as_quantum_mol(w.length, s.e.irrad))

sun.data <- as_tibble(sun.midday.data)
sun.spct <- as.source_spct(sun.data)
sun.spct <- trim_spct(sun.spct, low.limit = 280, fill = 0, use.hinges = FALSE)
setWhenMeasured(sun.spct, ymd_hms("2010-06-22 9:51:00"), tz = "UTC")
setWhereMeasured(sun.spct, geocode.kumpula)
setWhatMeasured(sun.spct, "sunlight, simulated")
comment(sun.spct) <- "Simulated solar spectrum based on real weather conditions.\nData author: Dr. Anders Lindfors\nFinnish Meteorological Institute.\nSee help file for references."

save(sun.spct, file = "data/sun.spct.rda")
save(sun.data, file = "data/sun.data.rda")
rm(sun.data, sun.spct, sun.midday.data)

sun.daily.data <- read.table("data-raw/sun_20120601_cum.hel.txt", col.names = c("w.length","s.e.irrad"))
sun.daily.data$w.length <- sun.daily.data$w.length / 10.0
sun.daily.data$s.q.irrad <- with(sun.daily.data, as_quantum_mol(w.length, s.e.irrad))
sun.daily.data <- as_tibble(sun.daily.data)
sun.daily.spct <- as.source_spct(sun.daily.data, time.unit = "day")
sun.daily.spct <- trim_spct(sun.daily.spct, low.limit = 280, fill = 0, use.hinges = FALSE)
setWhenMeasured(sun.daily.spct, ymd("2012-06-01"), tz = "UTC")
setWhereMeasured(sun.daily.spct, geocode.kumpula)
setWhatMeasured(sun.daily.spct, "sunlight, simulated")
comment(sun.daily.spct) <- "Total daily spectral exposure estimated from hourly simulations of the solar spectrum\nbased real weather conditions.\nData author: Dr. Anders Lindfors\nFinnish Meteorological Institute.\nSee help file for references."

save(sun.daily.spct, file = "data/sun.daily.spct.rda")
save(sun.daily.data, file = "data/sun.daily.data.rda")
rm(sun.daily.spct, sun.daily.data)
