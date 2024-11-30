.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Documentation at https://docs.r4photobiology.info/")
}

.onLoad <- function(libname, pkgname) {
  options(photobiology.verbose = getOption("verbose"))
}
