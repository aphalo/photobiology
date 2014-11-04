library(photobiology)
library(photobiologyWavebands)
library(microbenchmark)

data(sun.spct)
attach(sun.spct)

test.irradiance <- function(w.band=new_waveband(400,700)) {
  microbenchmark(irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=TRUE, use.cached.mult=FALSE),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=TRUE, use.cached.mult=TRUE),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=FALSE, use.cached.mult=TRUE),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=TRUE, use.cached.mult=FALSE, use.hinges=TRUE),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=FALSE, use.cached.mult=FALSE, use.hinges=TRUE),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.hinges=TRUE),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.hinges=FALSE)
)
}

test.irradiance()

test.irradiance(PAR())

test.irradiance(CIE())

test.irradiance(DNA.N())

test.irradiance(new_waveband(400,700))

test.irradiance(new_waveband(400,700,hinges=numeric(0)))

test.irradiance(new_waveband(400,700,hinges=NULL))

photon_irradiance(w.length, s.e.irrad, PAR(), use.hinges=TRUE)

photon_irradiance(w.length, s.e.irrad, PAR(), use.hinges=FALSE)

energy_irradiance(w.length, s.e.irrad, GEN.G(), use.hinges=TRUE)

energy_irradiance(w.length, s.e.irrad, GEN.G(), use.hinges=FALSE)

Rprof("irradiance.out")
for (i in 1:10000){
irradiance(w.length, s.e.irrad, new_waveband(400,700),"photon", check.spectrum=FALSE, use.cached.mult=TRUE)
}
Rprof(NULL)
summaryRprof("irradiance.out")

unlink("irradiance.out")

Rprof("irradiance.out")
for (i in 1:10000){
  irradiance(w.length, s.e.irrad, new_waveband(400,700, wb.name="my.test"),"photon", check.spectrum=FALSE, use.cached.mult=FALSE)
}
Rprof(NULL)
summaryRprof("irradiance.out")

unlink("irradiance.out")

Rprof("irradiance.out")
for (i in 1:10000){
  irradiance(w.length, s.e.irrad, new_waveband(400,700, hinges=numeric(0)),"photon", check.spectrum=FALSE, use.cached.mult=TRUE)
}
Rprof(NULL)
summaryRprof("irradiance.out")

unlink("irradiance.out")

Rprof("irradiance.out")
for (i in 1:10000){
  irradiance(w.length, s.e.irrad, DNA.N(), unit.out="energy", check.spectrum=TRUE, use.cached.mult=TRUE, use.hinges=TRUE)
}
Rprof(NULL)
summaryRprof("irradiance.out")

unlink("irradiance.out")

