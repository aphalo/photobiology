
# normalize ---------------------------------------------------------------

#' Normalize spectral data
#'
#' These functions return a spectral object of the same class as the one
#' supplied as argument but with the spectral data normalized to 1.o a certain
#' wavelength.
#'
#' @param x An R object
#' @param ... not used in current version
#' @return A new object of the same class as \code{x}.
#' @export normalize
#' @note Accepted values for \code{norm} vary depending on the class of
#'   \code{x}
#' @family rescaling functions
#'
normalize <- function(x, ...) UseMethod("normalize")

#' @describeIn normalize Default for generic function
#'
#' @export
normalize.default <- function(x, ...) {
  warning("'normalize' is not defined for objects of class ", class(x)[1])
  return(x)
}


#' @keywords internal
#'
normalize_spct <- function(spct, range, norm, var.name) {
  stopifnot(is.any_spct(spct), !is.null(var.name), length(var.name) == 1, var.name %in% names(spct))
  if (is.null(range) || is.na(range)) {
    range <- spct
  }
  tmp.spct <- trim_spct(spct, range)
  # rescaling needed
  if (!is.null(norm)) {
    if (is.character(norm)) {
      if (norm %in% c("max", "maximum")) {
        idx <- which.max(tmp.spct[[var.name]])
      } else if (norm %in% c("min", "minimum")) {
        idx <- which.min(tmp.spct[[var.name]])
      } else {
        warning("Invalid character '", norm, "'value in 'norm'")
        idx <- NA
      }
      scale.factor <- 1 / tmp.spct[idx, var.name, drop = TRUE]
      norm <- tmp.spct[idx, "w.length", drop = TRUE]
    } else if (is.numeric(norm)) {
      if (norm >= min(tmp.spct) && norm <= max(tmp.spct)) {
        tmp.spct <- tmp.spct[ , c("w.length", var.name)]
        class(tmp.spct) <- class(spct)
        scale.factor <- 1 / interpolate_spct(spct = tmp.spct, w.length.out = norm)[ , eval(var.name)]
      } else {
        warning("'norm = ", norm, "' value outside spectral data range of ",
                round(min(tmp.spct), 1), " to ", round(max(tmp.spct), 1), " (nm)")
        scale.factor <- NA
      }
    } else {
      stop("'norm' should be numeric or character")
    }
  } else {
    return(spct)
  }
  spct[[var.name]] <- spct[ , var.name, drop = TRUE] * scale.factor
  spct <- setNormalized(spct, norm)
  spct
}

#' @describeIn normalize Normalize a \code{source_spct} object.
#'
#' @param range An R object on which \code{range()} returns a numeric vector of
#'   length 2 with the limits of a range of wavelengths in nm, with min annd max
#'   wavelengths (nm)
#' @param norm numeric Normalization wavelength (nm) or character string "max",
#'   or "min" for normalization at the corresponding wavelngth, or "integral" or
#'   "mean" for rescaling by dividing by these values.
#' @param unit.out character Allowed values "energy", and "photon",
#'   or its alias "quantum"
#'
#' @export
#'
normalize.source_spct <- function(x,
                                  ...,
                                  range = NULL,
                                  norm = "max",
                                  unit.out = getOption("photobiology.radiation.unit", default="energy")) {
  if (unit.out == "energy") {
    return(normalize_spct(spct = q2e(x, action = "replace"),
                          range = range,
                          norm = norm,
                          var.name = "s.e.irrad"))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(normalize_spct(spct = e2q(x, action = "replace"),
                          range = range,
                          norm = norm,
                          var.name = "s.q.irrad"))
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' @describeIn normalize Normalize a response spectrum.
#'
#' @export
#'
normalize.response_spct <- function(x,
                                  ...,
                                  range = NULL,
                                  norm = "max",
                                  unit.out = getOption("photobiology.radiation.unit", default="energy")) {
  if (unit.out == "energy") {
    return(normalize_spct(spct = q2e(x, action = "replace"),
                          range = range,
                          norm = norm,
                          var.name = "s.e.response"))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(normalize_spct(spct = e2q(x, action = "replace"),
                          range = range,
                          norm = norm,
                          var.name = "s.q.response"))
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' @describeIn normalize Normalize a filter spectrum.
#'
#' @param qty.out character string  Allowed values are "transmittance", and
#'   "absorbance" indicating on which quantity to apply the normalization.
#'
#' @export
#'
normalize.filter_spct <- function(x,
                                  ...,
                                  range = NULL,
                                  norm = "max",
                                  qty.out = getOption("photobiology.filter.qty",
                                                      default="transmittance")) {
  if (qty.out == "transmittance") {
    return(normalize_spct(spct = A2T(x, action = "replace"),
                          range = range,
                          norm = norm,
                          var.name = "Tfr"))
  } else if (qty.out == "absorbance") {
    return(normalize_spct(spct = T2A(x, action = "replace"),
                          range = range,
                          norm = norm,
                          var.name = "A"))
  } else {
    stop("'qty.out ", qty.out, " is unknown")
  }
}

#' @describeIn normalize Normalize a reflector spectrum.
#'
#' @export
#'
normalize.reflector_spct <- function(x,
                                  ...,
                                  range = NULL,
                                  norm = "max",
                                  qty.out = NULL) {
    return(normalize_spct(spct = x,
                          range = range,
                          norm = norm,
                          var.name = "Rfr"))
}


# collections of spectra --------------------------------------------------


#' @describeIn normalize Normalize the members of a source_mspct object.
#'
#' @export
#'
normalize.source_mspct <- function(x,
                                   ...,
                                   range = NULL,
                                   norm = "max",
                                   unit.out = getOption("photobiology.radiation.unit",
                                                        default = "energy")) {
  msmsply(x,
          normalize,
          range = range,
          norm = norm,
          unit.out = unit.out,
          ...)

}

#' @describeIn normalize Normalize the members of a response_mspct object.
#'
#' @export
#'
normalize.response_mspct <- function(x,
                                     ...,
                                     range = NULL,
                                     norm = "max",
                                     unit.out = getOption("photobiology.radiation.unit",
                                                          default = "energy")) {
  msmsply(x,
          normalize,
          range = range,
          norm = norm,
          unit.out = unit.out,
          ...)

}

#' @describeIn normalize Normalize the members of a filter_mspct object.
#'
#' @export
#'
normalize.filter_mspct <- function(x,
                                   ...,
                                   range = NULL,
                                   norm = "max",
                                   qty.out = getOption("photobiology.filter.qty",
                                                       default = "transmittance")) {
  msmsply(x,
          normalize,
          range = range,
          norm = norm,
          qty.out = qty.out,
          ...)

}

#' @describeIn normalize Normalize the members of a reflector_mspct object.
#'
#' @export
#'
normalize.reflector_mspct <- function(x,
                                     ...,
                                     range = x,
                                     norm = "max",
                                     qty.out = NULL) {
  msmsply(x,
          normalize,
          range = range,
          norm = norm,
          qty.out = qty.out,
          ...)

}

# is_normalized function --------------------------------------------------

#' Query whether a generic spectrum has been normalized.
#'
#' This function tests a \code{generic_spct} object for an attribute that
#' signals whether the spectral data has been normalized or not after the object
#' was created.
#'
#' @param x An R object.
#'
#' @return A \code{logical} value. If \code{x} is not normalized or \code{x} is
#'   not a \code{generic_spct} object the value returned is \code{FALSE}.
#'
#' @export
#' @family rescaling functions
#'
is_normalized <- function(x) {
  if (!is.any_spct(x)) {
    return(NA)
  }
  spct.attr <- attr(x, "normalized", exact = TRUE)
  as.logical(!is.null(spct.attr) && as.logical(spct.attr))
}

# getNormalized -----------------------------------------------------------

#' Get the "normalized" attribute
#'
#' Funtion to read the "normalized" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#'
#' @return character or numeric or logical
#'
#' @note if x is not a \code{generic_spct} object, \code{NA} is returned
#'
#' @export
#' @family rescaling functions
#'
getNormalized <- function(x) {
  if (is.generic_spct(x)) {
    normalized <- attr(x, "normalized", exact = TRUE)
    if (is.null(normalized) || is.na(normalized)) {
      # need to handle objects created with old versions
      normalized <- FALSE
    }
    return(normalized[[1]])
  } else {
    return(NA)
  }
}

#' Set the "normalized" attribute
#'
#' Funtion to write the "normalized" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#' @param norm numeric or logical
#'
#' @note if x is not a \code{generic_spct} object, x is not modified.
#'
#' @export
#' @family rescaling functions
#'
setNormalized <- function(x, norm = FALSE) {
  name <- substitute(x)
  if (is.any_spct(x) && (is.na(norm) || norm)) {
    attr(x, "normalized") <- norm
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

