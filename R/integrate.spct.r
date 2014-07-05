##' Returns the integrals of the data variables in a spectrum.
##'
##' This function gives the result of integrating spectral data over
##' wavelengths.
##'
##' @usage integrate_spct(spct)
##'
##' @param spct an object of class "generic.spct"
##'
##' @return one or more numeric values with no change in scale factor: e.g. [W m-2 nm-1] -> [W m-2]
##' @keywords manip misc
##' @export
##' @examples
##' data(sun.data)
##' sun.spct <- setGenSpct(sun.data)
##' integrate_spct(sun.spct)
integrate_spct <- function(spct) {
  names.spct <- names(spct)
  names.data <-names.spct[names.spct != "w.length"]
  comment.spct <- comment(spct)
  first.iter <- TRUE
  integrals <- numeric(length(names.data))
  for (data.col in names.data) {
    integrals <- c(integrals, integrate_irradiance(spct[["w.length"]], spct[[eval(data.col)]]))
  }
  names(integrals) <- names.data
  setattr(integrals, "comment", comment.spct)
  return(integrals)
}
