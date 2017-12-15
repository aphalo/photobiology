library(photobiology)

clear.spct <- filter_spct(w.length = c(100, 101, 4999, 5000),
                          Tfr = rep(1, 4),
                          Tfr.type = "internal")
setWhatMeasured(clear.spct, "theoretical fully transparent object")
opaque.spct <- filter_spct(w.length = c(100, 101, 4999, 5000),
                           Tfr = rep(0, 4),
                           Tfr.type = "internal")
setWhatMeasured(opaque.spct, "theoretical fully opaque object")

library(photobiologyFilters)

polyester.spct <- filters.mspct$Autostat_CT5_125um
setWhatMeasured(polyester.spct, "clear polyester film, 125um thick")
yellow_gel.spct <- filters.mspct$Canary_Supergel312
setWhatMeasured(yellow_gel.spct, "yellow theatrical 'gel', Rosco supergel no. 312, 'canary yellow'")

save(clear.spct, file = "./data/clear.spct.rda")
save(opaque.spct, file = "./data/opaque.spct.rda")
save(polyester.spct, file = "./data/polyester.spct.rda")
save(yellow_gel.spct, file = "./data/yellow.gel.rda")
