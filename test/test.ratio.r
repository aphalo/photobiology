library(photobiology)
library(photobiologyUV)
library(photobiologyVIS)
library(microbenchmark)

data(sun.data)
attach(sun.data)

test.waveband_ratio <- function(w.band.num=new_waveband(315,400),w.band.denom=new_waveband(400,700)) {
  microbenchmark(waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=TRUE, use.cached.mult=FALSE),
                 waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=TRUE, use.cached.mult=TRUE),
                 waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=FALSE, use.cached.mult=TRUE),
                 waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.hinges=TRUE),
                 waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.hinges=FALSE),
                 waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=TRUE, use.cached.mult=FALSE, use.cpp.code=FALSE),
                 waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=TRUE, use.cached.mult=TRUE, use.cpp.code=FALSE),
                 waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.cpp.code=FALSE),
                 waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.cpp.code=FALSE, use.hinges=TRUE),
                 waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.cpp.code=FALSE, use.hinges=FALSE))
}

test.waveband_ratio()

test.waveband_ratio(UVB(),PAR())

test.waveband_ratio(Red("Smith"),Far_red())


wb.num <- UVB()
wb.denom <- PAR()

Rprof("ratio.out")
for (i in 1:1000){
waveband_ratio(w.length, s.e.irrad, wb.num, wb.denom,"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.cpp.code=TRUE)
}
Rprof(NULL)
summaryRprof("ratio.out")

unlink("ratio.out")