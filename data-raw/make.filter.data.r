library(photobiology)

clear.spct <- filter_spct(w.length = c(100, 101, 4999, 5000),
                          Tfr = rep(1, 4),
                          Tfr.type = "internal")
opaque.spct <- filter_spct(w.length = c(100, 101, 4999, 5000),
                           Tfr = rep(0, 4),
                           Tfr.type = "internal")

library(photobiologyFilters)

polyester.spct <- filters.mspct$Autostat_CT5_125um
yellow_gel.spct <- filters.mspct$Canary_Supergel312

save(clear.spct, file = "./data/clear.spct.rda")
save(opaque.spct, file = "./data/opaque.spct.rda")
save(polyester.spct, file = "./data/polyester.spct.rda")
save(yellow_gel.spct, file = "./data/yellow.gel.rda")
