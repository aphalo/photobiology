library(photobiology)
setwd("raw.data")
sun.midday.data <- read.table("sun_20100622_midday.txt", col.names=c("time_min","w.length","s.e.irrad"))
sun.midday.data$time_min <- NULL
sun.midday.data$w.length <- sun.midday.data$w.length / 10.0
sun.midday.data$s.e.irrad <- sun.midday.data$s.e.irrad / 1e3
sun.midday.data$s.q.irrad <- with(sun.midday.data, as_quantum_mol(w.length, s.e.irrad))

sun.spct <- sun.dt <- sun.data <- data.table(sun.midday.data)
setSourceSpct(sun.spct)
setDT(sun.dt)

setwd("../data")

save(sun.spct, sun.data, sun.dt, file="sun.spct.rda")
rm(sun.dt, sun.data, sun.spct, sun.midday.data)

setwd("..")

setwd("raw.data")

sun.daily.spct <- read.table("sun_20120601_cum.hel.txt", col.names=c("w.length","s.e.irrad"))
setSourceSpct(sun.daily.spct)
sun.daily.spct[ , w.length := w.length * 1e-1]
# sun.daily.spct[ , s.e.irrad := s.e.irrad * 1e-3]
e2q(sun.daily.spct)

sun.daily.data <- sun.daily.dt <- sun.daily.spct

setwd("../data")

save(sun.daily.spct, sun.daily.data, sun.daily.dt, file="sun.daily.spct.rda")
rm(sun.daily.spct, sun.daily.data, sun.daily.dt)

setwd("..")
