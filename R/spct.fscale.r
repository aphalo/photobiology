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
#' @param f character string "mean" or "total" for scaling so taht this summary
#'   value becomes 1 for the returned object, or the name of a function taking
#'   \code{x} as first argument and returning a numeric value.
#' @param unit.out character Alowed values "energy", and "photon", or its alias
#'   "quantum"
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

#' fscale a spectrum
#'
#' These functions return a spectral object of the same class as the one
#' supplied as argument but with the spectral data scaled.
#'
#' @param spct generic_spct The spectrum to be normalized
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
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
  tmp.spct <- tmp.spct[ , c("w.length", var.name)]
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
  out.spct <- spct
  out.spct[[var.name]] <- out.spct[[var.name]] / summary.value
  class(out.spct) <- class(spct)
  comment(out.spct) <- comment(spct)
  setScaled(out.spct, TRUE)
  setTimeUnit(out.spct, getTimeUnit(spct))
  setTfrType(out.spct, getTfrType(spct))
  out.spct
}

# is_scaled function ----------------------------------------------------

#' Query whether a generic spectrum has been scaled
#'
#' This function tests a \code{generic_spct} object for an attribute that
#' signals whether the spectral data has been rescled or not after the object
#' was created.
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

# getScaled -----------------------------------------------------------

#' Get the "scaled" attribute
#'
#' Funtion to read the "scaled" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#'
#' @return logical
#'
#' @note if x is not a \code{filter_spct} object, \code{NA} is returned
#'
#' @export
#' @family Rfr attribute functions
#'
getScaled <- function(x) {
  if (is.any_spct(x)) {
    scaled <- attr(x, "scaled", exact = TRUE)
    if (is.null(scaled) || is.na(scaled)) {
      # need to handle objects created with old versions
      scaled <- FALSE
    }
    return(scaled[[1]])
  } else {
    return(NA)
  }
}

#' Set the "scaled" attribute
#'
#' Funtion to write the "scaled" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#' @param scaled logical
#'
#' @note if x is not a \code{generic_spct} object, x is not modified.
#'   attribute set.
#'
#' @export
#' @family rescaling functions
#'
setScaled <- function(x, scaled = FALSE) {
  if (is.na(scaled) || scaled) {
    attr(x, "scaled") <- scaled
  }
  return(x)
}

