library(photobiology)
setwd("data-raw")
sun.midday.data <- read.table("sun_20100622_midday.txt", col.names=c("time_min","w.length","s.e.irrad"))
sun.midday.data$time_min <- NULL
sun.midday.data$w.length <- sun.midday.data$w.length / 10.0
sun.midday.data$s.e.irrad <- sun.midday.data$s.e.irrad / 1e3
sun.midday.data$s.q.irrad <- with(sun.midday.data, as_quantum_mol(w.length, s.e.irrad))

sun.spct <- sun.data <- sun.midday.data
setSourceSpct(sun.spct)

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
sun.daily.spct <- sun.daily.data

setSourceSpct(sun.daily.spct, time.unit="day")

setwd("../data")

save(sun.daily.spct, file="sun.daily.spct.rda")
save(sun.daily.data, file="sun.daily.data.rda")
rm(sun.daily.spct, sun.daily.data)

setwd("..")
