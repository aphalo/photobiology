#' Rescale a spectrum.
#'
#' This function returns a spectral object of the same class as the one supplied
#' as argument but with the spectral data rescaled to 1.0 at a .
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#' The object returned is a ggplot objects, and can be further manipulated.
#'
#' @usage rescale_spct(spct,
#'                     range,
#'                     var.name,
#'                     f,
#'                     ...)
#'
#' @param spct a generic.spct object
#' @param range an R object on which range() returns a vector of length 2,
#' with min annd max wavelengths (nm)
#' @param var.name the name of the variable to Rescale
#' @param f a summary function to be applied to \code{spct}
#' @param ... other arguments passed to f()
#'
#' @return a new object of the same class as \code{spct}.
#'
#' @keywords internal
#'

rescale_spct <- function(spct, range, var.name, f, ...) {
  stopifnot(is.any.spct(spct), !is.null(var.name), length(var.name) == 1, var.name %in% names(spct))
  tmp.spct <- trim_spct(spct, range)
  tmp.spct <- tmp.spct[ , .SD, .SDcols = c("w.length", var.name)]
  # rescaling needed
  if (!is.null(f)) {
    if (is.character(f)) {
      if (f %in% c("mean", "average")) {
        summary.value <- average_spct(tmp.spct)
      } else if (f %in% c("total", "integral")) {
        summary.value <- integrate_spct(tmp.spct)
      } else {
        warning("Invalid character '", f, "'value in 'f'")
        summary.value <- NA_real_
      }
    } else if (is.function(f)) {
      summary.value <- f(tmp.spct, ...)
    } else {
      stop("'f' should be a function name or character")
    }
  } else {
    summary.value <- 1 # implemented in this way to ensure that all returned
    # values folow the same copy/reference semantics
  }
  out.spct <- copy(spct)
  out.spct[ , var.name := out.spct[ , unlist(.SD), .SDcols = var.name] / summary.value, with = FALSE]
  setattr(out.spct, "class", class(spct))
  setattr(out.spct, "comment", comment(spct))
  setattr(out.spct, "spct", c(attr(spct, "spct", exact = TRUE), rescaled = TRUE, f = f))
  out.spct
}


#' Generic function
#'
#' Function that returns an object with rescaled data.
#'
#' @param x an R object
#' @param ... not used in current version
#' @export Rescale
Rescale <- function(x, ...) UseMethod("Rescale")

#' Default for generic function
#'
#' Function that returns an object with rescaled data.
#'
#' @param x an R object
#' @param ... not used in current version
#' @export Rescale.default
Rescale.default <- function(x, ...) {
  return(x)
}

#' Rescale a source spectrum.
#'
#' This function returns a spectral object of the same class as the one supplied
#' as argument but with the spectral data rescaled to 1.0 at a given wavelength.
#'
#' @usage Rescale(x,
#'                ...,
#'                range = range(x),
#'                f = "mean",
#'                unit.out = getOption("photobiology.radiation.unit", default="energy"))
#'
#' @param spct a source.spct object
#' @param ... not used in current version
#' @param range an R object on which range() returns a vector of length 2,
#' with min annd max wavelengths (nm)
#' @param f numeric normalization wavelength (nm) or character string "mean",
#' or "total" for normalization at the corresponding wavelength.
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#'
#' @return a new object of the same class as \code{x}.
#'
#' @export
#'
Rescale.source.spct <- function(x,
                                ...,
                                range = x,
                                f = "mean",
                                unit.out = getOption("photobiology.radiation.unit", default="energy")) {
  if (unit.out == "energy") {
    return(rescale_spct(spct = q2e(x, action = "replace"),
                        range = range(range),
                        f = f,
                        var.name = "s.e.irrad"))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(rescale_spct(spct = e2q(x, action = "replace"),
                        range = range(range),
                        f = f,
                        var.name = "s.q.irrad"))
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' Rescale a response spectrum.
#'
#' This function returns a spectral object of the same class as the one supplied
#' as argument but with the spectral data rescaled to 1.0 at a given wavelength.
#'
#' @usage Rescale(x,
#'                ...,
#'                range = range(x),
#'                f = "mean",
#'                unit.out = getOption("photobiology.radiation.unit", default="energy"))
#'
#' @param spct a response.spct object
#' @param ... not used in current version
#' @param range an R object on which range() returns a vector of length 2,
#' with min annd max wavelengths (nm)
#' @param f numeric normalization wavelength (nm) or character string "mean",
#' or "total" for normalization at the corresponding wavelngth.
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#'
#' @return a new object of the same class as \code{x}.
#'
#' @export
#'
Rescale.response.spct <- function(x,
                                  ...,
                                  range = x,
                                  f = "mean",
                                  unit.out = getOption("photobiology.radiation.unit", default="energy")) {
  if (unit.out == "energy") {
    return(rescale_spct(spct = q2e(x, action = "replace"),
                        range = range(range),
                        f = f,
                        var.name = "s.e.response",
                        ...))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(rescale_spct(spct = e2q(x, action = "replace"),
                        range = range(range),
                        f = f,
                        var.name = "s.q.response",
                        ...))
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' Rescale a filter spectrum.
#'
#' This function returns a spectral object of the same class as the one supplied
#' as argument but with the spectral data rescaled to 1.0 at a given wavelength.
#'
#' @usage Rescale(x,
#'                ...,
#'                range = range(x),
#'                f = "mean",
#'                qty.out = getOption("photobiology.filter.qty", default="transmittance"))
#'
#' @param spct a filter.spct object
#' @param ... not used in current version
#' @param range an R object on which range() returns a vector of length 2,
#' with min annd max wavelengths (nm)
#' @param f numeric normalization wavelength (nm) or character string "mean",
#' or "total" for normalization at the corresponding wavelngth.
#' @param qty.out character string with allowed values "transmittance", and "absorbance"
#'
#' @return a new object of the same class as \code{x}.
#'
#' @export
#'
Rescale.filter.spct <- function(x,
                                ...,
                                range = x,
                                f = "mean",
                                qty.out = getOption("photobiology.filter.qty", default="transmittance")) {
  if (qty.out == "transmittance") {
    return(rescale_spct(spct = A2T(x, action = "replace"),
                        range = range(range),
                        f = f,
                        var.name = "Tfr",
                        ...))
  } else if (qty.out == "absorbance") {
    return(rescale_spct(spct = T2A(x, action = "replace"),
                        range = range(range),
                        f = f,
                        var.name = "A",
                        ...))
  } else {
    stop("'qty.out ", unit.out, " is unknown")
  }
}

#' Rescale a reflector spectrum.
#'
#' This function returns a spectral object of the same class as the one supplied
#' as argument but with the spectral data rescaled to 1.0 at a given wavelength.
#'
#' @usage Rescale(x,
#'                ...,
#'                range = range(x),
#'                f = "mean",
#'                qty.out = NULL)
#'
#' @param spct a reflector.spct object
#' @param ... not used in current version
#' @param range an R object on which range() returns a vector of length 2,
#' with min annd max wavelengths (nm)
#' @param f numeric normalization wavelength (nm) or character string "mean",
#' or "total" for normalization at the corresponding wavelngth.
#' @param qty.out ignored
#'
#' @return a new object of the same class as \code{x}.
#'
#' @export
#'
Rescale.reflector.spct <- function(x,
                                   ...,
                                   range = x,
                                   f = "mean",
                                   qty.out = NULL) {
  return(rescale_spct(spct = x,
                      range = range(range),
                      f = f,
                      var.name = "Rfr",
                      ...))
}
