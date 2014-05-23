library(photobiology)
setwd("raw.data")
sun.midday.data <- read.table("sun_20100622_midday.txt", col.names=c("time_min","w.length","s.e.irrad"))
sun.midday.data$time_min <- NULL
sun.midday.data$w.length <- sun.midday.data$w.length/10.0
sun.midday.data$s.e.irrad <- sun.midday.data$s.e.irrad / 1e3
sun.midday.data$s.q.irrad <- with(sun.midday.data, as_quantum_mol(w.length, s.e.irrad))
setwd("../data")

sun.dt <- sun.midday.data
setSourceSpct(sun.dt)
sun.data <- sun.dt
save(sun.dt, sun.data, file="sun.dt.rda")
rm(sun.dt, sun.data, sun.midday.data)
setwd("..")
