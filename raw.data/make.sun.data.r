setwd("raw.data")
sun.midday.data <- read.table("sun_20100622_midday.txt", col.names=c("time_min","w.length","s.e.irrad"))
sun.midday.data$time_min <- NULL
sun.midday.data$w.length <- sun.midday.data$w.length/10.0
sun.midday.data$s.e.irrad <- sun.midday.data$s.e.irrad / 1e3
setwd("../data")
save(sun.midday.data, file="sun.data")
setwd("..")
