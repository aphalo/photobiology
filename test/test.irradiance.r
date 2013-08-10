library(photobiology)
library(photobiologyUV)
library(microbenchmark)

data(sun.data)
attach(sun.data)

test.irradiance <- function(w.band=new_waveband(400,700)) {
  microbenchmark(irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=TRUE, use.cached.mult=FALSE),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=TRUE, use.cached.mult=TRUE),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=FALSE, use.cached.mult=TRUE),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=TRUE, use.cached.mult=FALSE, use.cpp.code=FALSE),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=TRUE, use.cached.mult=TRUE, use.cpp.code=FALSE),
                 irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.cpp.code=FALSE))
}

test.irradiance()

test.irradiance(new_waveband(400,700))

test.irradiance(new_waveband(400,700,hinges=numeric(0)))

test.irradiance(new_waveband(400,700,hinges=NULL))

Rprof("irradiance.out")
irradiance(w.length, s.e.irrad, new_waveband(400,700),"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.cpp.code=FALSE)
Rprof(NULL)
summaryRprof("irradiance.out")

Rprof("irradiance.out")
test.irradiance()
Rprof(NULL)
summaryRprof("irradiance.out")

Rprof("irradiance.out")
for (i in 1:1000){
irradiance(w.length, s.e.irrad, new_waveband(400,700),"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.cpp.code=TRUE)
}
Rprof(NULL)
summaryRprof("irradiance.out")

Rprof("irradiance.out")
for (i in 1:1000){
  irradiance(w.length, s.e.irrad, new_waveband(400,700, hinges=numeric(0)),"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.cpp.code=TRUE)
}
Rprof(NULL)
summaryRprof("irradiance.out")
unlink("irradiance.out")