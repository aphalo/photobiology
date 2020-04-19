library(photobiology)
library(ooacquire)

load(file = "./data-raw/ooacquire/led_desk201.spct.Rda")

white_led.raw_spct <- led_desk201.raw_spct$light
white_led.cps_spct <-
  s_irrad_corrected(led_desk201.raw_spct,
                    correction.method = MAYP11278_ylianttila.mthd, return.cps = TRUE)
white_led.source_spct <-
  s_irrad_corrected(led_desk201.raw_spct,
                    correction.method = MAYP11278_ylianttila.mthd)

white_led.source_spct <- despike(white_led.source_spct)
white_led.source_spct <- fshift(white_led.source_spct)
white_led.source_spct <- smooth_spct(white_led.source_spct)
white_led.source_spct <- clean(white_led.source_spct)

trimInstrDesc(white_led.raw_spct, c("-", "w", "inst.calib"))
trimInstrDesc(white_led.cps_spct, c("-", "w", "inst.calib"))
trimInstrDesc(white_led.source_spct, c("-", "w", "inst.calib"))

getInstrDesc(white_led.raw_spct)
getInstrDesc(white_led.cps_spct)
getInstrDesc(white_led.source_spct)

getInstrSettings(white_led.raw_spct)
getInstrSettings(white_led.cps_spct)
getInstrSettings(white_led.source_spct)

save(white_led.raw_spct, white_led.cps_spct, white_led.source_spct,
     file = "./data/white-led-spct.rda")

