library(photobiology)

photodiode.spct <- read.csv("data-raw/photodiode/hamamatsu_G6262.csv", col.names = c("w.length","s.e.response"))
setResponseSpct(photodiode.spct, time.unit = "second")

photodiode.spct <- trim_wl(photodiode.spct, range = c(300, 580), fill = NULL)
comment(photodiode.spct) <- "Spectral response of GaAsP photodiode.\nHamamatsu G6262\nSource: GaAsP photodiodes datasheet,\n document at http://www.hamamatsu.com/, Hamamatsu (2011), Hamamatsu City, Japan.\nResponse expressed in A/W."
setWhatMeasured(photodiode.spct, "GaAsP photodiode")

ccd.spct <- read.csv("data-raw/photodiode/hamamatsu_S10420-1.csv", col.names = c("w.length","s.q.response"))
setResponseSpct(ccd.spct, time.unit = "second")

ccd.spct <- trim_wl(ccd.spct, range = c(200, 1100), fill = NULL)
comment(ccd.spct) <- "Spectral response of CCD.\nHamamatsu S10420-1\nSource: CCD image sensors datasheet,\n document at http://www.hamamatsu.com/, Hamamatsu (2014), Hamamatsu City, Japan.\nResponse expressed quantum efficiency (as fraction of one)."
setWhatMeasured(ccd.spct, "CCD linear image sensor")

save(photodiode.spct, file = "data/photodiode.spct.rda")
save(ccd.spct, file = "data/ccd.spct.rda")
rm(photodiode.spct, ccd.spct)
