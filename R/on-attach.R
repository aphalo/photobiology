.onAttach <- function(libname, pkgname) {
  packageStartupMessage("For news on the R for Photobiology packages, please, visit https://www.r4photobiology.info/")
}

.onLoad <- function(libname, pkgname) {
  options(photobiology.verbose = getOption("verbose"))
}
