library(photobiology)

clear.spct <- filter_spct(w.length = c(100, 101, 4999, 5000),
                          Tfr = rep(1, 4),
                          Tfr.type = "internal",
                          Rfr.constant = 0.00,
                          thickness = 1,
                          attenuation.mode = "absorption")
setWhatMeasured(clear.spct, "theoretical fully transparent object")

opaque.spct <- filter_spct(w.length = c(100, 101, 4999, 5000),
                           Tfr = rep(0, 4),
                           Tfr.type = "internal",
                           Rfr.constant = 0.00,
                           thickness = 1,
                           attenuation.mode = "absorption")
setWhatMeasured(opaque.spct, "theoretical fully opaque object")

library(photobiologyFilters)

polyester.spct <- filters.mspct$McDermit_PET_Autostat_CT5_125um
setTfrType(polyester.spct, "total")
setWhatMeasured(polyester.spct, "clear polyester film")
setFilterProperties(polyester.spct,
                    Rfr.constant = 0.07, # guess
                    thickness = 125e-6,
                    attenuation.mode = "absorption")

yellow_gel.spct <- filters.mspct$Rosco_Canary_Supergel_no312
setTfrType(polyester.spct, "total")
setWhatMeasured(yellow_gel.spct, "yellow theatrical 'gel', Rosco supergel no. 312, 'canary yellow'")
setFilterProperties(yellow_gel.spct,
                    Rfr.constant = 0.07, # guess
                    thickness = 85e-6,
                    attenuation.mode = "absorption")

save(clear.spct, file = "./data/clear.spct.rda")
save(opaque.spct, file = "./data/opaque.spct.rda")
save(polyester.spct, file = "./data/polyester.spct.rda")
save(yellow_gel.spct, file = "./data/yellow.gel.rda")
