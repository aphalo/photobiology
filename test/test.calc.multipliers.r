library(photobiologyUV)

test.calc_multipliers <- function(w.band=new_waveband(400,700)) {
  data(sun.data)
  ptm0 <- proc.time()
  for (i in 1:10000){
    with(sun.data, calc_multipliers(w.length, w.band,"photon", use.cached.mult=FALSE))  
  }
  ptm1 <- proc.time()
  for (i in 1:10000){
    with(sun.data, calc_multipliers(w.length, w.band,"photon", use.cached.mult=TRUE))  
  }
  ptm2 <- proc.time()
  print(ptm0)
  print(ptm1 - ptm0)
  print(ptm2 - ptm1)
  ptm0 <- proc.time()
  for (i in 1:10000){
    with(sun.data, calc_multipliers(w.length, w.band,"photon", use.cached.mult=FALSE))  
  }
  ptm1 <- proc.time()
  for (i in 1:10000){
    with(sun.data, irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=TRUE, use.cached.mult=FALSE))  
  }
  ptm1 <- proc.time()
  for (i in 1:10000){
    with(sun.data, irradiance(w.length, s.e.irrad, w.band,"photon", check.spectrum=FALSE, use.cached.mult=TRUE))  
  }
  ptm2 <- proc.time()
  print(ptm0)
  print(ptm1 - ptm0)
  print(ptm2 - ptm1)
}

test.calc_multipliers()
test.calc_multipliers(DNA.N())
test.calc_multipliers(CIE())
test.calc_multipliers(CIE(300))

