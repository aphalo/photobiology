.onAttach <- function(libname, pkgname) {
  packageStartupMessage("For news and documentation about the R for Photobiology packages, please, visit http://www.r4photobiology.info/")
}

.onLoad <- function(libname, pkgname) {
  options(photobiology.verbose = getOption("verbose"))
}
