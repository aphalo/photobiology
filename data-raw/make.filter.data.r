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

# we fetch the spectra without thinning
load("../photobiologyFilters/data-raw/rda/rosco.mspct.rda")
load("../photobiologyFilters/data-raw/rda/mcdermit.mspct.rda")

polyester.spct <- mcdermit.mspct$McDermit_PET_Autostat_CT5_125um

yellow_gel.spct <- rosco.mspct$Rosco_Canary_Supergel_no312

two_filters.mspct <- filter_mspct(list(polyester.spct, yellow_gel.spct))

two_filters.spct <- rbindspct(list(polyester.spct, yellow_gel.spct))

save(clear.spct, opaque.spct, polyester.spct, yellow_gel.spct,
     two_filters.spct, two_filters.mspct, file = "data/filter-data.rda")
