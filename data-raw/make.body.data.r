library(photobiology)

black_body.spct <- object_spct(w.length = c(100, 101, 4999, 5000),
                          Tfr = rep(0, 4), Rfr = rep(0, 4),
                          Tfr.type = "internal",
                          Rfr.type = "total")
white_body.spct <- object_spct(w.length = c(100, 101, 4999, 5000),
                               Tfr = rep(0, 4), Rfr = rep(1, 4),
                               Tfr.type = "internal",
                               Rfr.type = "total")
clear_body.spct <- object_spct(w.length = c(100, 101, 4999, 5000),
                               Tfr = rep(1, 4), Rfr = rep(0, 4),
                               Tfr.type = "internal",
                               Rfr.type = "total")

save(black_body.spct, file = "./data/black_body.spct.rda")
save(white_body.spct, file = "./data/white_body.spct.rda")
save(clear_body.spct, file = "./data/clear_body.spct.rda")
