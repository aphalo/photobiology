#' Normalize a spectrum.
#'
#' This function returns a spectral object of the same class as the one supplied
#' as argument but with the spectral data normalized to 1.0 at a .
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#' The object returned is a ggplot objects, and can be further manipulated.
#'
#' @usage normalize_spct(spct,
#'                       range,
#'                       norm,
#'                       var.name)
#'
#' @param spct a generic.spct object
#' @param range an R object on which range() returns a vector of length 2,
#' with min annd max wavelengths (nm)
#' @param norm numeric normalization wavelength (nm) or character string "max",
#' or "min" for normalization at the corresponding wavelngth, or "integral"
#' or "mean" for rescaling by dividing by these values.
#' @param var.name the name of the variable to normalize
#'
#' @return a new object of the same class as \code{spct}.
#'
#' @keywords internal
#'

normalize_spct <- function(spct, range, norm, var.name) {
  stopifnot(is.any.spct(spct), !is.null(var.name), length(var.name) == 1, var.name %in% names(spct))
  tmp.spct <- trim_spct(spct, range)
  # rescaling needed
  if (!is.null(norm)) {
    if (is.character(norm)) {
      if (norm %in% c("max", "maximum")) {
        idx <- which.max(tmp.spct[ , unlist(.SD), .SDcols = var.name])
      } else if (norm %in% c("min", "minimum")) {
        idx <- which.min(tmp.spct[ , unlist(.SD), .SDcols = var.name])
      } else {
        warning("Invalid character '", norm, "'value in 'norm'")
        idx <- NA
      }
      scale.factor <- 1 / as.numeric(tmp.spct[idx, .SD, .SDcols = var.name])
      norm <- tmp.spct[idx, w.length]
    } else if (is.numeric(norm)) {
      if (norm >= min(tmp.spct) && norm <= max(tmp.spct)) {
        tmp.spct <- tmp.spct[ , .SD, .SDcols = c("w.length", var.name)]
        setattr(tmp.spct, "class", class(spct))
        scale.factor <- 1 / interpolate_spct(tmp.spct, norm)[ , unlist(.SD), .SDcols = var.name]
      } else {
        warning("'norm = ", norm, "' value outside spectral data range of ",
                round(min(tmp.spct), 1), " to ", round(max(tmp.spct), 1), " (nm)")
        scale.factor <- NA
      }
    } else {
      stop("'norm' should be numeric or character")
    }
  } else {
    scale.factor <- 1 # implemented in this way to ensure that all returned
    # values folow the same copy/reference semantics
  }
  out.spct <- copy(spct)
  out.spct[ , var.name := out.spct[ , unlist(.SD), .SDcols = var.name] * scale.factor, with = FALSE]
  setattr(out.spct, "class", class(spct))
  setattr(out.spct, "comment", comment(spct))
  setattr(out.spct, "normalized", norm)
  out.spct
}

#' Generic function
#'
#' Function that returns an object with normalized data.
#'
#' @param x an R object
#' @param ... not used in current version
#' @export normalize
normalize <- function(x, ...) UseMethod("normalize")

#' Default for generic function
#'
#' Function that returns an object with normalized data.
#'
#' @param x an R object
#' @param ... not used in current version
#' @export normalize.default
normalize.default <- function(x, ...) {
  return(x / max(x, ...))
}

#' Normalize a source spectrum.
#'
#' This function returns a spectral object of the same class as the one supplied
#' as argument but with the spectral data normalized to 1.0 at a given wavelength.
#'
#' @usage normalize(x,
#'                  ...,
#'                  range = range(x),
#'                  norm = "max",
#'                  unit.out = getOption("photobiology.radiation.unit", default="energy"))
#'
#' @param x a source.spct object
#' @param ... not used in current version
#' @param range an R object on which range() returns a vector of length 2,
#' with min annd max wavelengths (nm)
#' @param norm numeric normalization wavelength (nm) or character string "max",
#' or "min" for normalization at the corresponding wavelength.
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#'
#' @return a new object of the same class as \code{x}.
#'
#' @export
#'
normalize.source.spct <- function(x,
                                  ...,
                                  range = x,
                                  norm = "max",
                                  unit.out = getOption("photobiology.radiation.unit", default="energy")) {
  if (unit.out == "energy") {
    return(normalize_spct(spct = q2e(x, action = "replace"),
                          range = range(range),
                          norm = norm,
                          var.name = "s.e.irrad"))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(normalize_spct(spct = e2q(x, action = "replace"),
                          range = range(range),
                          norm = norm,
                          var.name = "s.q.irrad"))
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' Normalize a response spectrum.
#'
#' This function returns a spectral object of the same class as the one supplied
#' as argument but with the spectral data normalized to 1.0 at a given wavelength.
#'
#' @usage normalize(x,
#'                  ...,
#'                  range = range(x),
#'                  norm = "max",
#'                  unit.out = getOption("photobiology.radiation.unit", default="energy"))
#'
#' @param x a response.spct object
#' @param ... not used in current version
#' @param range an R object on which range() returns a vector of length 2,
#' with min annd max wavelengths (nm)
#' @param norm numeric normalization wavelength (nm) or character string "max",
#' or "min" for normalization at the corresponding wavelngth.
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#'
#' @return a new object of the same class as \code{x}.
#'
#' @export
#'
normalize.response.spct <- function(x,
                                  ...,
                                  range = x,
                                  norm = "max",
                                  unit.out = getOption("photobiology.radiation.unit", default="energy")) {
  if (unit.out == "energy") {
    return(normalize_spct(spct = q2e(x, action = "replace"),
                          range = range(range),
                          norm = norm,
                          var.name = "s.e.response"))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(normalize_spct(spct = e2q(x, action = "replace"),
                          range = range(range),
                          norm = norm,
                          var.name = "s.q.response"))
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' Normalize a filter spectrum.
#'
#' This function returns a spectral object of the same class as the one supplied
#' as argument but with the spectral data normalized to 1.0 at a given wavelength.
#'
#' @usage normalize(x,
#'                  ...,
#'                  range = range(x),
#'                  norm = "max",
#'                  qty.out = getOption("photobiology.filter.qty", default="transmittance"))
#'
#' @param spct a filter.spct object
#' @param ... not used in current version
#' @param range an R object on which range() returns a vector of length 2,
#' with min annd max wavelengths (nm)
#' @param norm numeric normalization wavelength (nm) or character string "max",
#' or "min" for normalization at the corresponding wavelngth.
#' @param qty.out character string with allowed values "transmittance", and "absorbance"
#'
#' @return a new object of the same class as \code{x}.
#'
#' @export
#'
normalize.filter.spct <- function(x,
                                  ...,
                                  range = x,
                                  norm = "max",
                                  qty.out = getOption("photobiology.filter.qty", default="transmittance")) {
  if (qty.out == "transmittance") {
    return(normalize_spct(spct = A2T(x, action = "replace"),
                          range = range(range),
                          norm = norm,
                          var.name = "Tfr"))
  } else if (qty.out == "absorbance") {
    return(normalize_spct(spct = T2A(x, action = "replace"),
                          range = range(range),
                          norm = norm,
                          var.name = "A"))
  } else {
    stop("'qty.out ", unit.out, " is unknown")
  }
}

#' Normalize a reflector spectrum.
#'
#' This function returns a spectral object of the same class as the one supplied
#' as argument but with the spectral data normalized to 1.0 at a given wavelength.
#'
#' @usage normalize(x,
#'                  ...,
#'                  range = range(x),
#'                  norm = "max",
#'                  qty.out = NULL)
#'
#' @param x a reflector.spct object
#' @param ... not used in current version
#' @param range an R object on which range() returns a vector of length 2,
#' with min annd max wavelengths (nm)
#' @param norm numeric normalization wavelength (nm) or character string "max",
#' or "min" for normalization at the corresponding wavelngth.
#' @param qty.out ignored
#'
#' @return a new object of the same class as \code{x}.
#'
#' @export
#'
normalize.reflector.spct <- function(x,
                                  ...,
                                  range = x,
                                  norm = "max",
                                  qty.out = NULL) {
    return(normalize_spct(spct = x,
                          range = range(range),
                          norm = norm,
                          var.name = "Rfr"))
}

#' Query whether a generic spectrum has been normalized.
#'
#' This function returns TRUE if x is a generic.spct and it has been normalized
#' by means of function \code{normalize}.
#'
#' @usage is.rescaled(x)
#'
#' @param x a generic.spct object
#'
#' @export
#'
is.normalized <- function(x) {
  if (!is.any.spct(x)) {
    return(NA)
  }
  spct.attr <- attr(x, "normalized", exact = TRUE)
  as.logical(!is.null(spct.attr) && as.logical(spct.attr))
}

