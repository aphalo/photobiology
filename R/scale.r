# fscale methods ---------------------------------------------------------

#' Rescale a spectrum using a summary function
#'
#' These functions return a spectral object of the same class as the one supplied
#' as argument but with the spectral data rescaled.
#'
#' @param x An R object
#' @param ... additonal named arguments passed down to \code{f}.
#' @export fscale
#' @family rescaling functions
#'
fscale <- function(x, ...) UseMethod("fscale")

#' @describeIn fscale Default for generic function
#'
#' @export
#'
#' @return a new object of the same class as \code{x}.
#'
fscale.default <- function(x, ...) {
  warning("'fscale' is not defined for objects of class ", class(spct)[1])
  return(x)
}

#' @describeIn fscale
#'
#' @param range An R object on which \code{range()} returns a numeric vector of
#'   length 2 with the limits of a range of wavelengths in nm, with min annd max
#'   wavelengths (nm)
#' @param f numeric Normalization wavelength (nm) or character string "mean",
#' or "total" for normalization at the corresponding wavelength.
#' @param unit.out character Alowed values "energy", and "photon", or its alias "quantum"
#'
#' @export
#'
fscale.source_spct <- function(x,
                               range = x,
                               f = "mean",
                               unit.out = getOption("photobiology.radiation.unit", default="energy"),
                               ...) {
  if (unit.out == "energy") {
    return(fscale_spct(spct = q2e(x, action = "replace"),
                       range = range(range),
                       f = f,
                       var.name = "s.e.irrad"))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(fscale_spct(spct = e2q(x, action = "replace"),
                       range = range(range),
                       f = f,
                       var.name = "s.q.irrad"))
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' @describeIn fscale
#'
#' @export
#'
fscale.response_spct <- function(x,
                                 range = x,
                                 f = "mean",
                                 unit.out = getOption("photobiology.radiation.unit", default="energy"),
                                 ...) {
  if (unit.out == "energy") {
    return(fscale_spct(spct = q2e(x, action = "replace"),
                       range = range(range),
                       f = f,
                       var.name = "s.e.response",
                       ...))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(fscale_spct(spct = e2q(x, action = "replace"),
                       range = range(range),
                       f = f,
                       var.name = "s.q.response",
                       ...))
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' @describeIn fscale
#'
#' @param qty.out character Allowed values "transmittance", and "absorbance"
#'
#' @export
#'
fscale.filter_spct <- function(x,
                               range = x,
                               f = "mean",
                               qty.out = getOption("photobiology.filter.qty", default="transmittance"),
                               ...) {
  if (qty.out == "transmittance") {
    return(fscale_spct(spct = A2T(x, action = "replace"),
                       range = range(range),
                       f = f,
                       var.name = "Tfr",
                       ...))
  } else if (qty.out == "absorbance") {
    return(fscale_spct(spct = T2A(x, action = "replace"),
                       range = range(range),
                       f = f,
                       var.name = "A",
                       ...))
  } else {
    stop("'qty.out ", unit.out, " is unknown")
  }
}

#' @describeIn fscale
#'
#' @export
#'
fscale.reflector_spct <- function(x,
                                  range = x,
                                  f = "mean",
                                  qty.out = NULL,
                                  ...) {
  return(fscale_spct(spct = x,
                     range = range(range),
                     f = f,
                     var.name = "Rfr",
                     ...))
}

#' fscale a spectrum.
#'
#' These functions return a spectral object of the same class as the one supplied
#' as argument but with the spectral data scaled.
#'
#' @usage fscale_spct(spct,
#'                     range,
#'                     var.name,
#'                     f,
#'                     ...)
#'
#' @param spct generic_spct The spectrum to be normalized
#' @param range an R object on which range() returns a vector of length 2,
#' with min annd max wavelengths (nm)
#' @param var.name character The name of the variable to fscale
#' @param f function A summary function to be applied to \code{spct}
#' @param ... other arguments passed to f()
#'
#' @return a new object of the same class as \code{spct}.
#'
#' @keywords internal
#'
fscale_spct <- function(spct, range, var.name, f, ...) {
  stopifnot(is.any_spct(spct), !is.null(var.name), length(var.name) == 1, var.name %in% names(spct))
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
  setattr(out.spct, "scaled", TRUE)
  setTimeUnit(out.spct, getTimeUnit(spct))
  setTfrType(out.spct, getTfrType(spct))
  out.spct
}

# is_scaled function ----------------------------------------------------

#' Query whether a generic spectrum has been scaled.
#'
#' This function tests a \code{generic_spct} object for an attribute that
#' signals whether the spectral data has been rescled or not after the object
#' was created.
#'
#' @usage is_scaled(x)
#'
#' @param x An R object.
#'
#' @return A \code{logical} value. If \code{x} is not scaled or \code{x} is
#'   not a \code{generic_spct} object the value returned is \code{FALSE}.
#'
#' @export
#' @family rescaling functions
#'
is_scaled <- function(x) {
  if (!is.any_spct(x)) {
    return(NA)
  }
  spct.attr <- attr(x, "scaled", exact = TRUE)
  as.logical(!is.null(spct.attr) && as.logical(spct.attr))
}
