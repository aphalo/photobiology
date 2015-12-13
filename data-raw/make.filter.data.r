library(photobiology)

clear.spct <- filter_spct(w.length = c(100, 101, 4999, 5000),
                          Tfr = rep(1, 4))
opaque.spct <- filter_spct(w.length = c(100, 101, 4999, 5000),
                           Tfr = rep(0, 4))

library(photobiologyFilters)

polyester.spct <- mcdermit.mspct$Autostat_CT5_125um
yellow.gel.spct <- rosco.mspct$Canary_Supergel312

save(clear.spct, file = "./data/clear.spct.rda")
save(opaque.spct, file = "./data/opaque.spct.rda")
save(polyester.spct, file = "./data/polyester.spct.rda")
save(yellow.gel.spct, file = "./data/yellow.gel.rda")
