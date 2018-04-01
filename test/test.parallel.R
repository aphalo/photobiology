library(photobiology)
library(parallel)
library(foreach)
library(doParallel)
library(microbenchmark)

detectCores()

cluster <- makeCluster(4)
registerDoParallel(cluster)
# registerDoParallel(cores=4)

load("./test/test-data/spectra-mspct.rda")

many.mspct <- list(a = light.irrad_mspct,
                   b = light.irrad_mspct,
                   c = light.irrad_mspct,
                   d = light.irrad_mspct)

print(
  microbenchmark(
    foreach(i = 1:4) %do%
      q_irrad(many.mspct[[i]], .parallel = FALSE),
    times=5, unit="s")
)

print(
  microbenchmark(
    foreach(i = 1:4, .packages = "photobiology") %dopar%
      q_irrad(many.mspct[[i]], .parallel = FALSE),
    times=5, unit="s")
)
