# fscale methods ---------------------------------------------------------

#' Rescale a spectrum using a summary function
#'
#' These functions return a spectral object of the same class as the one supplied
#' as argument but with the spectral data rescaled.
#'
#' @param x An R object
#' @param ... additional named arguments passed down to \code{f}.
#'
#' @return A copy of \code{x} with the original spectral data values replaced
#'   with rescaled values, and the \code{"scaled"} attribute set to a list
#'   describing the scaling applied.
#'
#' @examples
#' fscale(sun.spct, f = "mean")
#' fscale(sun.spct, f = "mean", na.rm = TRUE)
#' fscale(sun.spct, f = sum)
#' fscale(sun.spct, f = function(x) {sum(x) / length(x)})
#'
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
  warning("'fscale()' is not defined for objects of class '", class(x)[1], "'.")
  return(x)
}

#' @describeIn fscale
#'
#' @param range An R object on which \code{range()} returns a numeric vector of
#'   length 2 with the limits of a range of wavelengths in nm, with min and max
#'   wavelengths (nm)
#' @param f character string "mean" or "total" for scaling so that this summary
#'   value becomes 1 for the returned object, or the name of a function taking
#'   \code{x} as first argument and returning a numeric value.
#' @param target numeric A constant used as target value for scaling.
#' @param unit.out character Allowed values "energy", and "photon", or its alias
#'   "quantum"
#'
#' @export
#'
fscale.source_spct <- function(x,
                               range = NULL,
                               f = "mean",
                               target = 1,
                               unit.out = getOption("photobiology.radiation.unit", default="energy"),
                               ...) {
  if (unit.out == "energy") {
    return(fscale_spct(spct = q2e(x, action = "replace"),
                       range = range,
                       f = f,
                       target = target,
                       col.names = "s.e.irrad",
                       ...))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(fscale_spct(spct = e2q(x, action = "replace"),
                       range = range,
                       f = f,
                       target = target,
                       col.names = "s.q.irrad",
                       ...))
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' @describeIn fscale
#'
#' @export
#'
fscale.response_spct <- function(x,
                                 range = NULL,
                                 f = "mean",
                                 target = 1,
                                 unit.out = getOption("photobiology.radiation.unit", default="energy"),
                                 ...) {
  if (unit.out == "energy") {
    return(fscale_spct(spct = q2e(x, action = "replace"),
                       range = range,
                       f = f,
                       target = target,
                       col.names = "s.e.response",
                       ...))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(fscale_spct(spct = e2q(x, action = "replace"),
                       range = range,
                       f = f,
                       target = target,
                       col.names = "s.q.response",
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
                               range = NULL,
                               f = "mean",
                               target = 1,
                               qty.out = getOption("photobiology.filter.qty",
                                                   default = "transmittance"),
                               ...) {
  if (qty.out == "transmittance") {
    return(fscale_spct(spct = A2T(x, action = "replace"),
                       range = range,
                       f = f,
                       target = target,
                       col.names = "Tfr",
                       ...))
  } else if (qty.out == "absorbance") {
    return(fscale_spct(spct = T2A(x, action = "replace"),
                       range = range,
                       f = f,
                       target = target,
                       col.names = "A",
                       ...))
  } else {
    stop("'qty.out ", qty.out, " is unknown")
  }
}

#' @describeIn fscale
#'
#' @export
#'
fscale.reflector_spct <- function(x,
                                  range = NULL,
                                  f = "mean",
                                  target = 1,
                                  qty.out = NULL,
                                  ...) {
  return(fscale_spct(spct = x,
                     range = range,
                     f = f,
                     target = target,
                     col.names = "Rfr",
                     ...))
}

#' @describeIn fscale
#'
#' @export
#'
fscale.raw_spct <- function(x,
                            range = NULL,
                            f = "mean",
                            target = 1,
                            ...) {
  return(fscale_spct(spct = x,
                     range = range,
                     f = f,
                     target = target,
                     col.names = grep("^counts", names(x), value = TRUE),
                     ...))
}

#' @describeIn fscale
#'
#' @export
#'
fscale.cps_spct <- function(x,
                            range = NULL,
                            f = "mean",
                            target = 1,
                            ...) {
  return(fscale_spct(spct = x,
                     range = range,
                     f = f,
                     target = target,
                     col.names = grep("^cps", names(x), value = TRUE),
                     ...))
}

#' @describeIn fscale
#'
#' @param col.names character vector containing the names of columns or
#'   variables to which to apply the scaling.
#'
#' @export
#'
fscale.generic_spct <- function(x,
                                range = NULL,
                                f = "mean",
                                target = 1,
                                col.names,
                                ...) {
  return(fscale_spct(spct = x,
                     range = range,
                     f = f,
                     target = target,
                     col.names = col.names,
                     ...))
}

# Collections of spectra --------------------------------------------------

#' @describeIn fscale
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
fscale.source_mspct <- function(x,
                                range = NULL,
                                f = "mean",
                                target = 1,
                                unit.out = getOption("photobiology.radiation.unit",
                                                     default = "energy"),
                                ...,
                                .parallel = FALSE,
                                .paropts = NULL) {
  msmsply(x,
          fscale,
          range = range,
          f = f,
          target = target,
          unit.out = unit.out,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn fscale
#'
#' @export
#'
fscale.response_mspct <- function(x,
                                  range = NULL,
                                  f = "mean",
                                  target = 1,
                                  unit.out = getOption("photobiology.radiation.unit",
                                                       default = "energy"),
                                  ...,
                                  .parallel = FALSE,
                                  .paropts = NULL) {
  msmsply(x,
          fscale,
          range = range,
          f = f,
          target = target,
          unit.out = unit.out,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn fscale
#'
#' @export
#'
fscale.filter_mspct <- function(x,
                                range = NULL,
                                f = "mean",
                                target = 1,
                                qty.out = getOption("photobiology.filter.qty",
                                                    default = "transmittance"),
                                ...,
                                .parallel = FALSE,
                                .paropts = NULL) {
  msmsply(x,
          fscale,
          range = range,
          f = f,
          target = target,
          qty.out = qty.out,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn fscale
#'
#' @export
#'
fscale.reflector_mspct <- function(x,
                                   range = NULL,
                                   f = "mean",
                                   target = 1,
                                   qty.out = NULL,
                                   ...,
                                   .parallel = FALSE,
                                   .paropts = NULL) {
  msmsply(x,
          fscale,
          range = range,
          f = f,
          target = target,
          qty.out = qty.out,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn fscale
#'
#' @export
#'
fscale.raw_mspct <- function(x,
                             range = NULL,
                             f = "mean",
                             target = 1,
                             ...,
                             .parallel = FALSE,
                             .paropts = NULL) {
  msmsply(x,
          fscale,
          range = range,
          f = f,
          target = target,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn fscale
#'
#' @export
#'
fscale.cps_mspct <- function(x,
                             range = NULL,
                             f = "mean",
                             target = 1,
                             ...,
                             .parallel = FALSE,
                             .paropts = NULL) {
  msmsply(x,
          fscale,
          range = range,
          f = f,
          target = target,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn fscale
#'
#' @export
#'
fscale.generic_mspct <- function(x,
                                 range = NULL,
                                 f = "mean",
                                 target = 1,
                                 col.names,
                                 ...,
                                 .parallel = FALSE,
                                 .paropts = NULL) {
  msmsply(x,
          fscale,
          range = range,
          f = f,
          target = target,
          col.names = col.names,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}


# PRIVATE -----------------------------------------------------------------

#' fscale a spectrum
#'
#' These functions return a spectral object of the same class as the one
#' supplied as argument but with the spectral data scaled.
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
fscale_spct <- function(spct, range, col.names, f, target, ...) {
  # re-scaling will wipe out any existing normalization
  if (is_normalized(spct)) {
    setNormalized(spct, norm = FALSE)
  }

  if (is.null(range) || all(is.na(range))) {
    range <- range(spct)
  } else {
    if (max(range) < min(spct)) {
      warning("'range' does not overlap spectral data, skipping scaling...")
      return(spct)
    }
    if (min(range) < min(spct)) {
      warning("'range' is only partly within spectral data, continuing scaling...")
    }
  }
  tmp.spct <- trim_spct(spct, range, byref = FALSE)
  multipliers <- numeric(length(col.names))
  i <- 0L
  for (col in col.names) {
    i <- i + 1L
    # rescaling needed
    if (!is.null(f)) {
      if (is.character(f)) {
        if (f %in% c("mean", "average")) {
          summary.value <- average_spct(tmp.spct[, c("w.length", col)])
        } else if (f %in% c("total", "integral")) {
          summary.value <- integrate_spct(tmp.spct[, c("w.length", col)])
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
      summary.value <- 1 # implemented in this way to ensure that all returned
      # values follow the same copy/reference semantics
      if (target != 1) {
        warning("No summary function supplied\n",
                "spectral values multiplied by 'target'")
      }
    }
    multipliers[i] <- target / summary.value
    spct[[col]] <- spct[[col]] * multipliers[i]
  }
  setScaled(spct, list(multiplier = multipliers, f = f, range = range, target = target))
  spct
}

# is_scaled function ----------------------------------------------------

#' Query whether a generic spectrum has been scaled
#'
#' This function tests a \code{generic_spct} object for an attribute that
#' signals whether the spectral data has been rescaled or not after the object
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
  if (!is.generic_spct(x) && !is.summary_generic_spct(x)) {
    return(NA)
  }
  spct.attr <- attr(x, "scaled", exact = TRUE)
  as.logical(!is.null(spct.attr) && as.logical(spct.attr[[1]]))
}

# getScaled -----------------------------------------------------------

#' Get the "scaled" attribute
#'
#' Function to read the "scaled" attribute of an existing generic_spct
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
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    scaled <- attr(x, "scaled", exact = TRUE)
    if (is.null(scaled) || (!is.list(scaled) && all(is.na(scaled)))) {
      # need to handle objects created with old versions
      scaled <- FALSE
    } else if (is.list(scaled)) {
      if (!"target" %in% names(scaled)) {
        # cater for objects scaled before version 0.9.12
        scaled <- c(scaled, list(target = 1))
      }
      if (!"range" %in% names(scaled)) {
        # cater for objects scaled before version 0.9.12
        scaled <- c(scaled, list(range = c(NA_real_, NA_real_)))
      }
    }
    return(scaled)
  } else {
    return(NA)
  }
}

#' Set the "scaled" attribute
#'
#' Function to write the "scaled" attribute of an existing generic_spct
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
  name <- substitute(x)
  if ((is.generic_spct(x) || is.summary_generic_spct(x)) && !is.null(scaled)) {
    attr(x, "scaled") <- scaled
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}
