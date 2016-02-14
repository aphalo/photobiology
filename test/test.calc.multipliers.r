library(photobiology)
library(photobiologyWavebands)
library(microbenchmark)

Sys.info()
Sys.time()
sessionInfo()

data(sun.data)
attach(sun.data)

test.calc_multipliers <- function(w.band=new_waveband(400,700)) {
microbenchmark(calc_multipliers(w.length, w.band,"photon", use.cached.mult = FALSE),
               calc_multipliers(w.length, w.band,"photon", use.cached.mult = TRUE))
}

test.calc_multipliers()
test.calc_multipliers(DNA.N())
test.calc_multipliers(CIE())
test.calc_multipliers(CIE(300))

test.irradiance <- function(w.band = new_waveband(400,700), unit.out = "energy") {
  microbenchmark(irradiance(w.length, s.e.irrad, w.band,unit.out,
                            check.spectrum = TRUE, use.cached.mult = FALSE),
                 irradiance(w.length, s.e.irrad, w.band,unit.out,
                            check.spectrum = TRUE, use.cached.mult = TRUE),
                 irradiance(w.length, s.e.irrad, w.band,unit.out,
                            check.spectrum = FALSE, use.cached.mult = TRUE),
                 irradiance(w.length, s.e.irrad, w.band,unit.out,
                            check.spectrum = TRUE, use.cached.mult = FALSE,
                            use.hinges = TRUE),
                 irradiance(w.length, s.e.irrad, w.band,unit.out,
                            check.spectrum = TRUE, use.cached.mult = TRUE,
                            use.hinges = TRUE),
                 irradiance(w.length, s.e.irrad, w.band,unit.out,
                            check.spectrum = FALSE, use.cached.mult = TRUE,
                            use.hinges = TRUE),
                 irradiance(w.length, s.e.irrad, w.band,unit.out,
                            check.spectrum = TRUE, use.cached.mult = FALSE,
                            use.hinges = FALSE),
                 irradiance(w.length, s.e.irrad, w.band,unit.out,
                            check.spectrum = TRUE, use.cached.mult = TRUE,
                            use.hinges = FALSE),
                 irradiance(w.length, s.e.irrad, w.band,unit.out,
                            check.spectrum = FALSE, use.cached.mult = TRUE,
                            use.hinges = FALSE))
}

test.irradiance(unit.out = "energy")
test.irradiance(DNA.N(), unit.out = "energy")
test.irradiance(CIE(), unit.out = "energy")
test.irradiance(CIE(300), unit.out = "energy")

test.irradiance(unit.out = "photon")
test.irradiance(DNA.N(), unit.out = "photon")
test.irradiance(CIE(), unit.out = "photon")
test.irradiance(CIE(300), unit.out = "photon")


test.integrate_irradiance <- function() {
  microbenchmark(integrate_xy(w.length, s.e.irrad))
}

test.integrate_irradiance()

detach(sun.data)

test.irrad <- function(w.band=new_waveband(400,700)) {
  microbenchmark(e_irrad(sun.spct, w.band, use.cached.mult = FALSE, use.hinges = NULL),
                 e_irrad(sun.spct, w.band, use.cached.mult = TRUE, use.hinges = NULL),                   e_irrad(sun.spct, w.band, use.cached.mult = FALSE, use.hinges = TRUE),
                 e_irrad(sun.spct, w.band, use.cached.mult = TRUE, use.hinges = TRUE),
                 e_irrad(sun.spct, w.band, use.cached.mult = TRUE, use.hinges = FALSE),
                 e_irrad(sun.spct, w.band, use.cached.mult = FALSE, use.hinges = FALSE))
}

test.irrad()
test.irrad(DNA.N())
test.irrad(CIE())
test.irrad(CIE(300))
test.irrad(UV_bands())
test.irrad(VIS_bands())
test.irrad(GEN_G(300))
test.irrad(DNA.N())

sun.spct <- q2e(sun.spct, action = "replace")

test.irrad()
test.irrad(DNA.N())
test.irrad(CIE())
test.irrad(CIE(300))
test.irrad(UV_bands())
test.irrad(VIS_bands())
test.irrad(GEN_G(300))
test.irrad(DNA.N())
