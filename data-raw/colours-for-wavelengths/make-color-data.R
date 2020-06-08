library(photobiology)

wls <- (100:4000)
CMF.color <- color_of(x = wls, type = "CMF")
CC.color <- color_of(x = wls, type = "CC")

wl_colors.spct <- generic_spct(w.length = wls,
                                        CMF = CMF.color,
                                        CC = CC.color)

save(wl_colors.spct, file = "./R/sysdata.rda")


