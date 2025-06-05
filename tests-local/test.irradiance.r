library(photobiology)
library(photobiologyWavebands)
library(microbenchmark)

data(sun.spct)

test.irrad <- function(w.band=new_waveband(400,700)) {
  microbenchmark(irrad(sun.spct, w.band),
                 irrad(sun.spct, w.band,"photon"),
                 irrad(sun.spct, w.band,"photon", use.cached.mult=FALSE),
                 irrad(sun.spct, w.band,"photon", use.cached.mult=TRUE),
                 irrad(sun.spct, w.band,"photon", use.cached.mult=FALSE,
                       use.hinges=TRUE),
                 irrad(sun.spct, w.band,"photon", use.cached.mult=TRUE,
                       use.hinges=TRUE),
                 irrad(sun.spct, w.band,"photon", use.cached.mult=TRUE,
                       use.hinges=FALSE)
  )
}

test.irrad()

test.irrad(CIE())

test.irrad(DNA.N())

dna <- DNA.N()

test.irrad(dna)

attach(sun.spct)

test.irradiance <- function(w.band=new_waveband(400,700)) {
  microbenchmark(irradiance(w.length, s.e.irrad, w.band,"photon"),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=TRUE, use.cached.mult=FALSE),
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

### 2015

Rprof("irradiance.out")
for (i in 1:5000){
e_irrad(sun.spct, CIE(), use.hinges = FALSE, wb.trim = TRUE, use.cached.mult = FALSE)
}
Rprof(NULL)
summaryRprof("irradiance.out")

Rprof("irradiance.out")
for (i in 1:5000){
  e_irrad(sun.spct, CIE())
}
Rprof(NULL)
summaryRprof("irradiance.out")


microbenchmark(q_irrad(sun.spct, Plant_bands()))

test.mspct <- source_mspct(list(A = sun.spct, B = sun.daily.spct, C = sun.daily.spct, D = sun.daily.spct))

microbenchmark(q_irrad(test.mspct, Plant_bands(),
                       use.hinges = FALSE, wb.trim = TRUE, use.cached.mult = TRUE))

microbenchmark(q_irrad(test.mspct, Plant_bands(),
                       use.hinges = TRUE, wb.trim = TRUE, use.cached.mult = TRUE))

microbenchmark(q_irrad(sun.spct, Plant_bands(),
                       use.hinges = FALSE, wb.trim = TRUE, use.cached.mult = TRUE))

microbenchmark(q_irrad(sun.spct, Plant_bands(),
                       use.hinges = TRUE, wb.trim = TRUE, use.cached.mult = TRUE))

microbenchmark(photon_irradiance(sun.spct$w.length, sun.spct$s.e.irrad, PAR(),
                                 use.hinges = FALSE, use.cached.mult = TRUE))

microbenchmark(photon_irradiance(sun.spct$w.length, sun.spct$s.e.irrad, PAR(),
                                 use.hinges = TRUE, use.cached.mult = TRUE))

options(photobiology.use.hinges = TRUE)

microbenchmark(photon_irradiance(sun.spct$w.length, sun.spct$s.e.irrad, PAR()))
microbenchmark(photon_irradiance(sun.spct$w.length, sun.spct$s.e.irrad, CIE()))
microbenchmark(photon_irradiance(sun.spct$w.length, sun.spct$s.e.irrad, Plant_bands()))
microbenchmark(q_irrad(sun.spct, PAR()))
microbenchmark(q_irrad(sun.spct, CIE()))
microbenchmark(q_irrad(sun.spct, Plant_bands()))
microbenchmark(q_irrad(test.mspct, PAR()))
microbenchmark(q_irrad(test.mspct, CIE()))
microbenchmark(q_irrad(test.mspct, Plant_bands()))

options(photobiology.use.hinges = FALSE)

microbenchmark(photon_irradiance(sun.spct$w.length, sun.spct$s.e.irrad, PAR()))
microbenchmark(photon_irradiance(sun.spct$w.length, sun.spct$s.e.irrad, CIE()))
microbenchmark(photon_irradiance(sun.spct$w.length, sun.spct$s.e.irrad, Plant_bands()))
microbenchmark(q_irrad(sun.spct, PAR()))
microbenchmark(q_irrad(sun.spct, CIE()))
microbenchmark(q_irrad(sun.spct, Plant_bands()))
microbenchmark(q_irrad(test.mspct, PAR()))
microbenchmark(q_irrad(test.mspct, CIE()))
microbenchmark(q_irrad(test.mspct, Plant_bands()))

