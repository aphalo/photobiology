library(photobiology)

setwd("data-raw/photodiode")
photodiode.spct <- read.csv("hamamatsu_G6262.csv", col.names = c("w.length","s.e.response"))
setResponseSpct(photodiode.spct, time.unit = "second")

photodiode.spct <- trim_wl(photodiode.spct, range = c(200, 800), fill = 0)
comment(photodiode.spct) <- "Spectral response of GaAsP photodiode.\nHamamatsu G6262\nSource: GaAsP photodiodes datasheet,\n document at http://www.hamamatsu.com/, Hamamatsu (2011), Hamamatsu City, Japan.\nResponse expressed in A/W."

ccd.spct <- read.csv("hamamatsu_S10420-1.csv", col.names = c("w.length","s.q.response"))
setResponseSpct(ccd.spct, time.unit = "second")

ccd.spct <- trim_wl(ccd.spct, range = c(200, 1200), fill = 0)
comment(ccd.spct) <- "Spectral response of CCD.\nHamamatsu S10420-1\nSource: CCD image sensors datasheet,\n document at http://www.hamamatsu.com/, Hamamatsu (2010), Hamamatsu City, Japan.\nResponse expressed quantum efficiency (as fraction of one)."

setwd("../../data")

save(photodiode.spct, file = "photodiode.spct.rda")
save(ccd.spct, file = "ccd.spct.rda")
rm(photodiode.spct, ccd.spct)

setwd("..")

