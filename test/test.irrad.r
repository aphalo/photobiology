library(photobiology)
library(photobiologyWavebands)
library(microbenchmark)

Sys.info()
Sys.time()
sessionInfo()

UV_bands <- UV_bands()
VIS_bands <- VIS_bands()
PAR <- PAR()
CIE <- CIE()
DNA.N <- DNA_N()
wb.sun <- waveband(sun.spct)
wb.50 <- waveband(c(400,450))
wb.200 <- waveband(c(400,600))
wb.400 <- waveband(c(400,800))

test.irrad.all <- function(w.band=wb.sun) {
  microbenchmark(e_irrad(sun.spct, w.band),
                 e_irrad(sun.spct, w.band,
                         use.cached.mult = FALSE,
                         use.hinges = FALSE,
                         wb.trim = FALSE,
                         idx = FALSE) )
}

test.irrad.all()
test.irrad.all(wb.sun)
test.irrad.all(wb.50)
test.irrad.all(wb.200)
test.irrad.all(wb.400)

test.irrad.all(DNA.N)

test.irrad.cache <- function(w.band=wb.sun) {
  microbenchmark(e_irrad(sun.spct, w.band),
                 e_irrad(sun.spct, w.band,
                         use.cached.mult = FALSE,
                         use.hinges = FALSE,
                         wb.trim = FALSE,
                         idx = FALSE),
                 e_irrad(sun.spct, w.band,
                         use.cached.mult = TRUE,
                         use.hinges = FALSE,
                         wb.trim = FALSE,
                         idx = FALSE),
                 e_irrad(sun.spct, w.band,
                         use.hinges = FALSE,
                         wb.trim = FALSE,
                         idx = FALSE))
}

test.irrad.cache()
test.irrad.cache(VIS_bands)
test.irrad.cache(CIE)
test.irrad.cache(DNA.N)

test.irrad.hinges <- function(w.band=wb.sun) {
  microbenchmark(e_irrad(sun.spct, w.band),
                 e_irrad(sun.spct, w.band,
                         use.cached.mult = FALSE,
                         use.hinges = NULL,
                         wb.trim = FALSE,
                         idx = FALSE),
                 e_irrad(sun.spct, w.band,
                         use.cached.mult = FALSE,
                         use.hinges = FALSE,
                         wb.trim = FALSE,
                         idx = FALSE),
                 e_irrad(sun.spct, w.band,
                         use.cached.mult = FALSE,
                         use.hinges = TRUE,
                         wb.trim = FALSE,
                         idx = FALSE),
                 e_irrad(sun.spct, w.band,
                         use.cached.mult = FALSE,
                         wb.trim = FALSE,
                         idx = FALSE))
}

test.irrad.hinges()
test.irrad.hinges(VIS_bands)
test.irrad.hinges(CIE)
test.irrad.hinges(DNA.N)

test.irrad.trim <- function(w.band=wb.sun) {
  microbenchmark(e_irrad(sun.spct, w.band),
                 e_irrad(sun.spct, w.band,
                         use.cached.mult = FALSE,
                         use.hinges = FALSE,
                         wb.trim = FALSE,
                         idx = FALSE),
                 e_irrad(sun.spct, w.band,
                         use.cached.mult = FALSE,
                         use.hinges = FALSE,
                         wb.trim = TRUE,
                         idx = FALSE),
                 e_irrad(sun.spct, w.band,
                         use.cached.mult = FALSE,
                         use.hinges = FALSE,
                         idx = FALSE))
}

test.irrad.trim()
test.irrad.trim(VIS_bands)
test.irrad.trim(CIE)
test.irrad.trim(DNA.N)

test.irrad.idx <- function(w.band=wb.sun) {
  microbenchmark(e_irrad(sun.spct, w.band),
                 e_irrad(sun.spct, w.band,
                         use.cached.mult = FALSE,
                         use.hinges = FALSE,
                         wb.trim = FALSE,
                         idx = FALSE),
                 e_irrad(sun.spct, w.band,
                         use.cached.mult = FALSE,
                         use.hinges = FALSE,
                         wb.trim = FALSE,
                         idx = TRUE),
                 e_irrad(sun.spct, w.band,
                         use.cached.mult = FALSE,
                         use.hinges = FALSE,
                         wb.trim = FALSE))
}


test.irrad.idx()
test.irrad.idx(VIS_bands)
test.irrad.idx(CIE)
test.irrad.idx(DNA.N)


