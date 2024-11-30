library(photobiology)
library(dplyr)
library(lubridate)
library(ggspectra)

geocode.kumpula <- data.frame(lat = 60.20911,
                              lon = 24.96474,
                              address = "Kumpula, Helsinki, FI",
                              stringsAsFactors = FALSE)

# midday spectrum ---------------------------------------------------------

sun.midday.data <- read.table("data-raw/sun_20100622_midday.txt", col.names = c("time_min","w.length","s.e.irrad"))
sun.midday.data$time_min <- NULL
sun.midday.data$w.length <- sun.midday.data$w.length / 10.0
sun.midday.data$s.e.irrad <- sun.midday.data$s.e.irrad / 1e3
sun.midday.data$s.q.irrad <- with(sun.midday.data, as_quantum_mol(w.length, s.e.irrad))

sun.data <- as.data.frame(sun.midday.data)
sun.spct <- as.source_spct(sun.data)
sun.spct <- trim_spct(sun.spct, low.limit = 280, fill = 0, use.hinges = FALSE)
setWhenMeasured(sun.spct, ymd_hms("2010-06-22 9:51:00"), tz = "UTC")
setWhereMeasured(sun.spct, geocode.kumpula)
setWhatMeasured(sun.spct, "sunlight, simulated")
setHowMeasured(sun.spct, "Simulated using 'libRadtran'. Average for one hour, centred approximately on local solar noon.")
comment(sun.spct) <- "Simulated solar spectrum obtained with 'libRadtran' radiation transfer model using satellite-based estimates of ozone column depth and ground-based measurements of cloud and aerosol depths for real local weather conditions.\nData author: Dr. Anders Lindfors\nFinnish Meteorological Institute.\nSee help file for references."

autoplot(sun.spct)

# daily spectrum ----------------------------------------------------------

sun.daily.data <- read.table("data-raw/sun_20120601_cum.hel.txt", col.names = c("w.length","s.e.irrad"))
sun.daily.data$w.length <- sun.daily.data$w.length / 10.0
sun.daily.data$s.q.irrad <- with(sun.daily.data, as_quantum_mol(w.length, s.e.irrad))
sun.daily.data <- as_tibble(sun.daily.data)
sun.daily.spct <- as.source_spct(sun.daily.data, time.unit = "day")
sun.daily.spct <- trim_spct(sun.daily.spct, low.limit = 280, fill = 0, use.hinges = FALSE)
setWhenMeasured(sun.daily.spct, ymd("2012-06-01"), tz = "UTC")
setWhereMeasured(sun.daily.spct, geocode.kumpula)
setWhatMeasured(sun.daily.spct, "sunlight, simulated")
setHowMeasured(sun.daily.spct, "Simulated using 'libRadtran'. Day-long spectrum computed by integration of simulated one hour averages through the whole day.")
comment(sun.daily.spct) <- "Simulated solar spectrum obtained with 'libRadtran' radiation transfer model using satellite-based estimates of ozone column depth and ground-based measurements of cloud and aerosol depths for real local weather conditions.\nData author: Dr. Anders Lindfors\nFinnish Meteorological Institute.\nSee help file for references."

sun_daily.spct <- sun.daily.spct
sun_daily.data <- sun.daily.data

autoplot(sun.daily.spct)

# evening time series -----------------------------------------------------

geocode.viiki <- data.frame(lat = 60.227,
                            lon = 24.018,
                            address = "Viikki, Helsinki, FI",
                            stringsAsFactors = FALSE)

load("data-raw/ooacquire/cosine.hour.9.spct.Rda")

# ensure we remove dependency on 'ooacquire'
cosine.hour.9.spct <- trimInstrDesc(cosine.hour.9.spct)
cosine.hour.9.spct <- trimInstrSettings(cosine.hour.9.spct)

where_measured(cosine.hour.9.spct) <- geocode.viiki
where_measured(cosine.hour.9.spct)

autoplot(cosine.hour.9.spct, facets = 3)

cosine.hour.9.spct <- despike(cosine.hour.9.spct)

cosine.hour.9.spct <- smooth_spct(cosine.hour.9.spct)

cosine.hour.9.spct <- clean(cosine.hour.9.spct)

cosine.hour.9.spct <- trim_wl(cosine.hour.9.spct, range = c(290, 1000))

sun_evening.spct <- subset(cosine.hour.9.spct, spct.idx %in% paste("time", 1:5, sep = ".0"))

autoplot(sun_evening.spct, facets = 3)

autoplot(sun_evening.spct, facets = 3, unit.out = "photon")

sun_evening.mspct <- subset2mspct(sun_evening.spct)

autoplot(sun_evening.mspct, facets = 3)

autoplot(sun_evening.mspct, facets = 3, unit.out = "photon")

# save data objects to .rda file ------------------------------------------

save(sun.daily.spct, sun.daily.data,
     sun_daily.spct, sun_daily.data,
     sun.spct, sun.data,
     sun_evening.spct, sun_evening.mspct,
     file = "data/sun.data.rda")

rm(list = ls(pattern = "*"))
