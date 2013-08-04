library(photobiologyUV)
library(microbenchmark)

data(sun.data)
attach(sun.data)

test.calc_multipliers <- function(w.band=new_waveband(400,700)) {
microbenchmark(calc_multipliers(w.length, w.band,"photon", use.cached.mult=FALSE),
               calc_multipliers(w.length, w.band,"photon", use.cached.mult=TRUE))
}

test.calc_multipliers()
test.calc_multipliers(DNA.N())
test.calc_multipliers(CIE())
test.calc_multipliers(CIE(300))

test.irradiance <- function(w.band=new_waveband(400,700)) {
  microbenchmark(irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=TRUE, use.cached.mult=FALSE),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=TRUE, use.cached.mult=TRUE),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=FALSE, use.cached.mult=TRUE))
}

test.irradiance()
test.irradiance(DNA.N())
test.irradiance(CIE())
test.irradiance(CIE(300))

test.integrate_irradiance <- function() {
  microbenchmark(integrate_irradiance(w.length, s.e.irrad),
                 integrateirr(w.length, s.e.irrad))
}

test.integrate_irradiance()

detach(sun.data)

