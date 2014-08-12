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
##' @usage interpolate_spct(spct, w.length.out=NULL, fill.value = NA, length.out=NULL)
##'
##' @param spct an object of class "generic.spct"
##' @param w.length.out numeric array of wavelengths (nm)
##' @param fill.value a value to be assigned to out of range wavelengths
##' @param length.out numeric value
##'
##' @details If \code{length.out} it is a numeric value, then gives the number of rows in the
##' output, if it is \code{NULL}, the values in the numeric vector \code{w.length.out} are used.
##' If both are not \code{NULL} then the range of \code{w.length.out} and \code{length.out} are
##' used to generate a vector of wavelength. A value of \code{NULL} for \code{fill} prevents
##' extrapolation.
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
##' interpolate_spct(sun.spct, c(400,500), length.out=201)
##'
interpolate_spct <- function(spct, w.length.out=NULL, fill.value=NA, length.out=NULL) {
  if (!is.null(length.out) && (is.na(length.out) || length.out < 1L) ) {
    return(spct[NA])
  }
  if (!is.null(length.out)  && length.out == 1L) {
    if (is.null(w.length.out)) {
      w.length.out <- midpoint(spct)
    } else {
      w.length.out <- midpoint(w.length.out)
    }
  }
  if (!is.null(length.out) && length.out > 1) {
    if (is.null(w.length.out) || length(w.length.out) < 2L) {
      w.length.out <- seq(min(spct), max(spct), length.out=length.out)
    } else {
      w.length.out <- seq(min(w.length.out), max(w.length.out), length.out=length.out)
    }
  } else if (is.null(w.length.out)) {
    # nothing to do
    return(spct)
  }
  names.spct <- names(spct)
  names.data <- names.spct[names.spct != "w.length"]
  comment.spct <- comment(spct)
  setkey(spct, w.length)
  max.spct <- max(spct)
  min.spct <- min(spct)
  max.wl.out <- max(w.length.out)
  min.wl.out <- min(w.length.out)
  if (min.spct > min.wl.out && min.spct < max.wl.out) w.length.out <- c(min.spct, w.length.out)
  if (max.spct < max.wl.out && max.spct > min.wl.out) w.length.out <- c(w.length.out, max.spct)
  if (is.null(fill.value)) {
    w.length.out <- w.length.out[w.length.out >= min.spct & w.length.out <= max.spct]
    if (length(w.length.out) == 0) {
      return(spct[NA])
    }
  }
  w.length.out <- unique(sort(w.length.out))
  new.spct <- data.table(w.length = w.length.out)

  for (data.col in names.data) {
    temp.values <-  with(spct, get(data.col))
    if (is.numeric(temp.values)) {
      new.values <- interpolate_spectrum(spct$w.length,
                                         temp.values,
                                         w.length.out,
                                         fill.value)
#      new.spct[ , as.character(eval(expression(data.col))) := new.values]
     new.spct[ , eval(data.col) := new.values]
    }
  }
  setattr(new.spct, "comment", comment(spct))
  setattr(new.spct, "class", class(spct))
  setkey(new.spct, w.length)
  return(new.spct)
}
