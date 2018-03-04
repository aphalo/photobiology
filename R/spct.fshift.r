# fshift methods ---------------------------------------------------------

#' Shift the scale of a spectrum using a summary function
#'
#' These functions return a spectral object of the same class as the one
#' supplied as argument but with the spectral data on a shifted scale. A range
#' of wavelengths is taken a reference (zero or another numeric constant) and a
#' summary is calculated for this waveband. The difference between the computed
#' and reference value are used to shift the scale so that these two values
#' match in the returned object.
#'
#' @param x An R object
#' @param ... additional named arguments passed down to \code{f}.
#'
#' @return A copy of \code{x} with the spectral data values replaced with values
#'   zero-shifted.
#'
#' @export fshift
#' @family rescaling functions
#'
fshift <- function(x, ...) UseMethod("fshift")

#' @describeIn fscale Default for generic function
#'
#' @export
#'
#' @return a new object of the same class as \code{x}.
#'
fshift.default <- function(x, ...) {
  warning("'fshift()' is not defined for objects of class '", class(x)[1], "'.")
  return(x)
}

#' @describeIn fshift
#'
#' @param range An R object on which \code{range()} returns a numeric vector of
#'   length 2 with the limits of a range of wavelengths in nm, with min and max
#'   wavelengths (nm)
#' @param f character string "mean", "min" or "max" for scaling so that this
#'   summary value becomes the origin of the spectral data scale in the returned
#'   object, or the name of a function taking \code{x} as first argument and
#'   returning a numeric value.
#' @param unit.out character Allowed values "energy", and "photon", or its alias
#'   "quantum"
#'
#' @export
#'
fshift.source_spct <-
  function(x,
           range = c(min(x), min(x) + 10),
           f = "mean",
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...) {
    if (unit.out == "energy") {
      return(fshift_spct(
        spct = q2e(x, action = "replace"),
        range = range,
        f = f,
        col.names = "s.e.irrad"
      ))
    } else if (unit.out %in% c("photon", "quantum")) {
      return(fshift_spct(
        spct = e2q(x, action = "replace"),
        range = range,
        f = f,
        col.names = "s.q.irrad"
      ))
    } else {
      stop("'unit.out ", unit.out, " is unknown")
    }
  }

#' @describeIn fshift
#'
#' @export
#'
fshift.response_spct <-
  function(x,
           range = c(min(x), min(x) + 10),
           f = "mean",
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...) {
    if (unit.out == "energy") {
      return(fshift_spct(
        spct = q2e(x, action = "replace"),
        range = range,
        f = f,
        col.names = "s.e.response",
        ...
      ))
    } else if (unit.out %in% c("photon", "quantum")) {
      return(fshift_spct(
        spct = e2q(x, action = "replace"),
        range = range,
        f = f,
        col.names = "s.q.response",
        ...
      ))
    } else {
      stop("'unit.out ", unit.out, " is unknown")
    }
  }

#' @describeIn fshift
#'
#' @param qty.out character Allowed values "transmittance", and "absorbance"
#'
#' @export
#'
fshift.filter_spct <- function(x,
                               range = c(min(x), min(x) + 10),
                               f = "min",
                               qty.out = getOption("photobiology.filter.qty",
                                                   default = "transmittance"),
                               ...) {
  if (qty.out == "transmittance") {
    return(fshift_spct(spct = A2T(x, action = "replace"),
                       range = range,
                       f = f,
                       col.names = "Tfr",
                       ...))
  } else if (qty.out == "absorbance") {
    return(fshift_spct(spct = T2A(x, action = "replace"),
                       range = range,
                       f = f,
                       col.names = "A",
                       ...))
  } else {
    stop("'qty.out ", qty.out, " is unknown")
  }
}

#' @describeIn fshift
#'
#' @export
#'
fshift.reflector_spct <- function(x,
                                  range = c(min(x), min(x) + 10),
                                  f = "min",
                                  qty.out = NULL,
                                  ...) {
  return(fshift_spct(spct = x,
                     range = range,
                     f = f,
                     col.names = "Rfr",
                     ...))
}

#' @describeIn fshift
#'
#' @export
#'
fshift.source_mspct <-
  function(x,
           range =  c(min(x), min(x) + 10),
           f = "mean",
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...) {
    msmsply(x,
            fshift,
            range = range,
            f = f,
            unit.out = unit.out,
            ...)
  }

#' @describeIn fshift
#'
#' @export
#'
fshift.raw_spct <- function(x,
                            range = c(min(x), min(x) + 10),
                            f = "mean",
                            qty.out = NULL,
                            ...) {
  return(fshift_spct(spct = x,
                     range = range,
                     f = f,
                     col.names = grep("^counts", names(x), value = TRUE),
                     ...))
}

#' @describeIn fshift
#'
#' @export
#'
fshift.cps_spct <- function(x,
                            range = c(min(x), min(x) + 10),
                            f = "mean",
                            qty.out = NULL,
                            ...) {
  return(fshift_spct(spct = x,
                     range = range,
                     f = f,
                     col.names = grep("^cps", names(x), value = TRUE),
                     ...))
}

#' @describeIn fshift
#'
#' @param col.names character vector containing the names of columns or
#'   variables to which to apply the scale shift.
#'
#' @export
#'
fshift.generic_spct <- function(x,
                                range = c(min(x), min(x) + 10),
                                f = "mean",
                                col.names,
                                ...) {
  return(fshift_spct(spct = x,
                     range = range,
                     f = f,
                     col.names = col.names,
                     ...))
}

# Collections of spectra --------------------------------------------------

#' @describeIn fshift
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
fshift.response_mspct <-
  function(x,
           range = c(min(x), min(x) + 10),
           f = "mean",
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            fshift,
            range = range,
            f = f,
            unit.out = unit.out,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @describeIn fshift
#'
#' @export
#'
fshift.filter_mspct <-
  function(x,
           range = c(min(x), min(x) + 10),
           f = "min",
           qty.out = getOption("photobiology.filter.qty",
                               default = "transmittance"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            fshift,
            range = range,
            f = f,
            qty.out = qty.out,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @describeIn fshift
#'
#' @export
#'
fshift.reflector_mspct <-
  function(x,
           range = c(min(x), min(x) + 10),
           f = "min",
           qty.out = NULL,
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            fshift,
            range = range,
            f = f,
            qty.out = qty.out,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @describeIn fshift
#'
#' @export
#'
fshift.raw_mspct <-
  function(x,
           range = c(min(x), min(x) + 10),
           f = "min",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            fshift,
            range = range,
            f = f,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @describeIn fshift
#'
#' @export
#'
fshift.cps_mspct <-
  function(x,
           range = c(min(x), min(x) + 10),
           f = "min",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            fshift,
            range = range,
            f = f,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @describeIn fshift
#'
#' @export
#'
fshift.generic_mspct <-
  function(x,
           range = c(min(x), min(x) + 10),
           f = "min",
           col.names,
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            fshift,
            range = range,
            f = f,
            col.names = col.names,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }


# PRIVATE -----------------------------------------------------------------

#' fshift a spectrum
#'
#' These function returns a spectral object of the same class as the one
#' supplied as argument but with the spectral data on a shifted scale.
#'
#' @param spct generic_spct The spectrum to be normalized
#' @param range an R object on which range() returns a vector of length 2, with
#'   min and max wavelengths (nm)
#' @param col.names character The name of the variable to fscale
#' @param f function A summary function to be applied to \code{spct}
#' @param ... other arguments passed to f()
#'
#' @return a new object of the same class as \code{spct}.
#'
#' @keywords internal
#'
fshift_spct <- function(spct, range, col.names, f, ...) {
  if (is.null(range) ||
      (!is.null(range) && max(range) < min(spct)) ||
      (!is.null(range) && min(range) < min(spct)) ) {
      warning("'range' does not fully overlap spectral data or is NULL, skipping fshifting...")
      return(spct)
  }
  tmp.spct <- trim_spct(spct, range, byref = FALSE)
  for (col in col.names) {
    # shifting needed
    if (!is.null(f)) {
      if (is.character(f)) {
        if (f %in% c("mean", "average")) {
          summary.value <- average_spct(tmp.spct[, c("w.length", col)])
        } else if (f %in% c("min", "minimum")) {
          summary.value <- min(tmp.spct[[col]])
        } else if (f %in% c("max", "maximum")) {
          summary.value <- max(tmp.spct[[col]])
        } else {
          warning("Invalid character '", f, "'value in 'f'")
          summary.value <- NA_real_
        }
      } else if (is.function(f)) {
        summary.value <- f(tmp.spct[, c("w.length", col)], ...)
        f <- "a user supplied R function"
      } else {
        stop("'f' should be a function name or character")
      }
    } else {
      summary.value <- 0
      # implemented in this way to ensure that all returned
      # values follow the same copy/reference semantics
    }
    spct[[col]] <- spct[[col]] - summary.value
  }
  spct
}
