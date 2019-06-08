.onAttach <- function(libname, pkgname) {
  packageStartupMessage("News at https://www.r4photobiology.info/")
}

.onLoad <- function(libname, pkgname) {
  options(photobiology.verbose = getOption("verbose"))
}
