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

##' Returns the data columns recalculated in a spectrum at a new set of wavelengths
##'
##' This function gives the result of interpolating spectral data from one set of
##' wavelengths to a new one.
##'
##' @usage interpolate_spct(spct, w.length.out, fill.value = NA)
##'
##' @param spct an object of class "generic.spct"
##' @param w.length.out numeric array of wavelengths (nm)
##' @param fill.value a value to be assigned to out of range wavelengths
##'
##' @note The default \code{fill.value = NA} fills extrpolated values with NA. Giving NULL as
##' argument for \code{fill.value} deletes wavelengths outside the input data range from the
##' returned spectrum. A numerical value can be also be provided as fill. This function calls
##' \code{interpolate_spectrum} for each non-wavelength column in the input spectra object.
##'
##' @return one or more numeric values with no change in scale factor: e.g. [W m-2 nm-1] -> [W m-2]
##' @keywords manip misc
##' @export
##' @examples
##' data(sun.spct)
##' interpolate_spct(sun.spct, 400:500, NA)
##' interpolate_spct(sun.spct, 400:500, NULL)
##' interpolate_spct(sun.spct, seq(200, 1000, by=0.1), 0)
##'
interpolate_spct <- function(spct, w.length.out, fill.value = NA) {
  setkey(spct, w.length)
  names.spct <- names(spct)
  names.data <-names.spct[names.spct != "w.length"]
  comment.spct <- comment(spct)
  first.iter <- TRUE
  if (is.null(fill.value)) {
    w.length.out <- w.length.out[w.length.out >= min(spct$w.length) & w.length.out <= max(spct$w.length)]
    w.length.out <- (c(min(spct$w.length), w.length.out, max(spct$w.length)))
  }
  w.length.out <- unique(sort(w.length.out))
  new.spct <- data.table(w.length = w.length.out)

  for (data.col in names.data) {
    temp.values <-  with(spct, get(data.col))
    new.values <- interpolate_spectrum(spct$w.length,
                                       temp.values,
                                       w.length.out,
                                       fill.value)
   new.spct[ , eval(expression(data.col)) := new.values]
  }
  setattr(new.spct, "comment", comment(spct))
  setattr(new.spct, "class", class(spct))
  setkey(new.spct, w.length)
  return(new.spct)
}
