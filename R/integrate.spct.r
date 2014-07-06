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
##' data(sun.spct)
##' integrate_spct(sun.spct)
integrate_spct <- function(spct) {
  names.spct <- names(spct)
  names.data <-names.spct[names.spct != "w.length"]
  comment.spct <- comment(spct)
  first.iter <- TRUE
  integrals <- NULL
  for (data.col in names.data) {
    integrals <- c(integrals, integrate_irradiance(spct[["w.length"]], spct[[eval(data.col)]]))
  }
  names(integrals) <- gsub("^s.", x = names.data, replacement = "")
  setattr(integrals, "comment", comment.spct)
  return(integrals)
}
##' Returns the average of the data variables in a spectrum.
##'
##' This function gives the result of integrating spectral data over
##' wavelengths and dividing by the spread or span of wavelengths.
##'
##' @usage average_spct(spct)
##'
##' @param spct an object of class "generic.spct"
##'
##' @return one or more numeric values with no change in scale factor: e.g. [W m-2 nm-1] -> [W m-2]
##' @keywords manip misc
##' @export
##' @examples
##' data(sun.spct)
##' average_spct(sun.spct)
average_spct <- function(spct) {
  return(integrate_spct(spct) / (max(spct) - min(spct)))
}
