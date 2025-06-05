library(photobiology)
library(photobiologyWavebands)
library(microbenchmark)

data(sun.spct)
attach(sun.spct)

test.waveband_ratio <- function(w.band.num=new_waveband(315,400),w.band.denom=new_waveband(400,700)) {
  microbenchmark(waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=TRUE, use.cached.mult=FALSE),
                 waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=TRUE, use.cached.mult=TRUE),
                 waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=FALSE, use.cached.mult=TRUE),
                 waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.hinges=TRUE),
                 waveband_ratio(w.length, s.e.irrad, w.band.num,w.band.denom,"photon", check.spectrum=FALSE, use.cached.mult=TRUE, use.hinges=FALSE)
  )
}

test.waveband_ratio()

test.waveband_ratio(UVB(),PAR())

test.waveband_ratio(Red("Smith"),Far_red("Smith"))


wb.num <- UVB()
wb.denom <- PAR()

Rprof("ratio.out")
for (i in 1:1000){
waveband_ratio(w.length, s.e.irrad, wb.num, wb.denom,"photon", check.spectrum=FALSE, use.cached.mult=TRUE)
}
Rprof(NULL)
summaryRprof("ratio.out")

unlink("ratio.out")
