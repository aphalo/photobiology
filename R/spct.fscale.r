# fscale methods ---------------------------------------------------------

#' Rescale a spectrum using a summary function
#'
#' These methods return a spectral object of the same class as the one
#' supplied as argument but with the spectral data rescaled based on a summary
#' function \code{f} applied over a specific \code{range} of wavelengths and a
#' \code{target} value for the summary value. When the object contains
#' multiple spectra, the rescaling is applied separately to each spectrum.
#'
#' @param x An R object
#' @param ... additional named arguments passed down to \code{f}.
#'
#' @details After scaling, calling the function passed as argument to \code{f}
#'   with the scaled spectrum as argument, will return the value passed as
#'   argument to
#'   \code{target}. \strong{The default for \code{set.scaled} depends
#'   dynamically on the value passed to \code{target}.} Sometimes we rescale a
#'   spectrum to a "theoretical" value for the summary, while in other cases we
#'   rescale the spectrum to a real-world target value of, e.g., a reference
#'   energy irradiance. In the first case we say that the data are expressed in
#'   relative units, while in the second case we retain actual physical units.
#'   To indicate this, the default argument for `set.scaled` is TRUE when
#'   \code{target == 1}, assuming the first of these two situations, and false
#'   otherwise, assuming the second situation. These defaults can be overriden
#'   with an explicit \code{logical} argument passed to \code{set.scaled}.
#'   Scaling overrides any previous normalization with the spectrum tagged as
#'   not normalized.
#'
#'   Method \code{fscale} is implemented for \code{solute_spct} objects but
#'   as the spectral data stored in them are a description of an intensive
#'   property of a substance, scaling is unlikely to useful. To represent
#'   solutions of specific concentrations of solutes, \code{filter_spct} objects
#'   should be used instead.
#'
#' @section Important changes: Metadata describing the rescaling operation are
#'   stored in an attribute only if \code{set.scaled = TRUE} is passed to the
#'   call. The exact format and data stored in the attribute \code{"scaled"} has
#'   changed during the development history of the package. Spectra re-scaled
#'   with earlier versions will lack some information. To obtain the metadata in
#'   a consistent format irrespective of this variation use accessor
#'   \code{getScaling()}, which fills missing fields with \code{NA}.
#'
#' @return A copy of the object passed as argument to \code{x} with the original
#'   spectral data values replaced with rescaled values, and the \code{"scaled"}
#'   attribute set to a list describing the scaling applied.
#'
#' @examples
#'
#' fscale(sun.spct)
#' fscale(sun.spct, f = "mean") # same as default
#' fscale(sun.spct, f = "mean", na.rm = TRUE)
#' fscale(sun.spct, range = c(400, 700)) # default is whole spectrum
#' fscale(sun.spct, f = "e_irrad", range = c(400, 700))
#' s400.spct <- fscale(sun.spct,
#'                     f = e_irrad,
#'                     range = c(400, 700),
#'                     target = 400) # a target in W m-2
#' s400.spct
#' e_irrad(s400.spct, c(400, 700))
#'
#' @export
#'
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
  x
}

#' @describeIn fscale
#'
#' @param range numeric. An R object on which \code{range()} returns a numeric
#'   vector of length 2 with the limits of a range of wavelengths in nm, with
#'   min and max wavelengths (nm)
#' @param f character string. "mean" or "total" for scaling so that this summary
#'   value becomes 1 for the returned object, or the name of a function taking
#'   \code{x} as first argument and returning a numeric value.
#' @param target numeric A constant used as target value for scaling.
#' @param unit.out character. Allowed values "energy", and "photon", or its alias
#'   "quantum".
#' @param set.scaled logical or NULL Flag indicating if the data is to be marked
#'   as "scaled" or not.
#'
#' @export
#'
fscale.source_spct <- function(x,
                               range = NULL,
                               f = "mean",
                               target = 1,
                               unit.out = getOption("photobiology.radiation.unit", default="energy"),
                               set.scaled = target == 1,
                               ...) {
  if (getMultipleWl(x) > 1L) {
    # brute force and slow approach, unsuitable for long time series
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      fscale.source_mspct(x = mspct,
                          range = range,
                          f = f,
                          target = target,
                          unit.out = unit.out,
                          set.scaled = set.scaled,
                          ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  if (unit.out == "energy") {
    fscale_spct(spct = q2e(x, action = "replace"),
                range = range,
                f = f,
                target = target,
                set.scaled = set.scaled,
                col.names = "s.e.irrad",
                ...)
  } else if (unit.out %in% c("photon", "quantum") ) {
    fscale_spct(spct = e2q(x, action = "replace"),
                range = range,
                f = f,
                target = target,
                set.scaled = set.scaled,
                col.names = "s.q.irrad",
                ...)
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
                                 set.scaled = target == 1,
                                 ...) {
  if (getMultipleWl(x) > 1L) {
    # brute force and slow approach, unsuitable for long time series
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      fscale.response_mspct(x = mspct,
                            range = range,
                            f = f,
                            target = target,
                            unit.out = unit.out,
                            set.scaled = set.scaled,
                            ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  if (unit.out == "energy") {
    fscale_spct(spct = q2e(x, action = "replace"),
                range = range,
                f = f,
                target = target,
                set.scaled = set.scaled,
                col.names = "s.e.response",
                ...)
  } else if (unit.out %in% c("photon", "quantum") ) {
    fscale_spct(spct = e2q(x, action = "replace"),
                range = range,
                f = f,
                target = target,
                set.scaled = set.scaled,
                col.names = "s.q.response",
                ...)
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' @describeIn fscale
#'
#' @param qty.out character. Allowed values "transmittance", and "absorbance".
#'
#' @export
#'
fscale.filter_spct <- function(x,
                               range = NULL,
                               f = "mean",
                               target = 1,
                               qty.out = getOption("photobiology.filter.qty",
                                                   default = "transmittance"),
                               set.scaled = target == 1,
                               ...) {
  if (getMultipleWl(x) > 1L) {
    # brute force and slow approach, unsuitable for long time series
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      fscale.filter_mspct(x = mspct,
                          range = range,
                          f = f,
                          target = target,
                          qty.out = qty.out,
                          set.scaled = set.scaled,
                          ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  if (qty.out == "transmittance") {
    fscale_spct(spct = A2T(x, action = "replace"),
                range = range,
                f = f,
                target = target,
                set.scaled = set.scaled,
                col.names = "Tfr",
                ...)
  } else if (qty.out == "absorbance") {
    fscale_spct(spct = T2A(x, action = "replace"),
                range = range,
                f = f,
                target = target,
                set.scaled = set.scaled,
                col.names = "A",
                ...)
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
                                  set.scaled = target == 1,
                                  ...) {
  if (getMultipleWl(x) > 1L) {
    # brute force and slow approach, unsuitable for long time series
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      fscale.reflector_mspct(x = mspct,
                             range = range,
                             f = f,
                             target = target,
                             qty.out = qty.out,
                             set.scaled = set.scaled,
                             ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  fscale_spct(spct = x,
              range = range,
              f = f,
              target = target,
              col.names = "Rfr",
              set.scaled = set.scaled,
              ...)
}

#' @describeIn fscale
#'
#' @export
#'
fscale.solute_spct <- function(x,
                               range = NULL,
                               f = "mean",
                               target = 1,
                               qty.out = NULL,
                               set.scaled = target == 1,
                               ...) {
  if (getMultipleWl(x) > 1L) {
    # brute force and slow approach, unsuitable for long time series
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      fscale.solute_mspct(x = mspct,
                          range = range,
                          f = f,
                          target = target,
                          qty.out = qty.out,
                          set.scaled = set.scaled,
                          ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  col.name <- intersect(c("K.mole", "K.mass"), names(x))
  fscale_spct(spct = x,
              range = range,
              f = f,
              target = target,
              col.names = col.name,
              set.scaled = set.scaled,
              ...)
}

#' @describeIn fscale
#'
#' @export
#'
fscale.raw_spct <- function(x,
                            range = NULL,
                            f = "mean",
                            target = 1,
                            set.scaled = target == 1,
                            ...) {
  if (getMultipleWl(x) > 1L) {
    # brute force and slow approach, unsuitable for long time series
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      fscale.raw_mspct(x = mspct,
                       range = range,
                       f = f,
                       target = target,
                       set.scaled = set.scaled,
                       ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  fscale_spct(spct = x,
              range = range,
              f = f,
              target = target,
              set.scaled = set.scaled,
             col.names = grep("^counts", names(x), value = TRUE),
               ...)
}

#' @describeIn fscale
#'
#' @export
#'
fscale.cps_spct <- function(x,
                            range = NULL,
                            f = "mean",
                            target = 1,
                            set.scaled = target == 1,
                            ...) {
  if (getMultipleWl(x) > 1L) {
    # brute force and slow approach, unsuitable for long time series
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      fscale.cps_mspct(x = mspct,
                       range = range,
                       f = f,
                       target = target,
                       set.scaled = set.scaled,
                       ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  fscale_spct(spct = x,
              range = range,
              f = f,
              target = target,
              set.scaled = set.scaled,
              col.names = grep("^cps", names(x), value = TRUE),
              ...)
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
                                set.scaled = target == 1,
                                col.names,
                                ...) {
  if (getMultipleWl(x) > 1L) {
    # brute force and slow approach, unsuitable for long time series
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      fscale.generic_mspct(x = mspct,
                           range = range,
                           f = f,
                           target = target,
                           set.scaled = set.scaled,
                           col.names = col.names,
                           ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  fscale_spct(spct = x,
              range = range,
              f = f,
              target = target,
              set.scaled = set.scaled,
              col.names = col.names,
              ...)
}

# Collections of spectra --------------------------------------------------

#' @describeIn fscale
#'
#' @param .parallel	logical if TRUE, apply function in parallel, using parallel
#'   backend provided by foreach.
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
                                set.scaled = target == 1,
                                ...,
                                .parallel = FALSE,
                                .paropts = NULL) {

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          fscale.source_spct,
          range = range,
          f = f,
          target = target,
          set.scaled = set.scaled,
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
                                  set.scaled = target == 1,
                                  ...,
                                  .parallel = FALSE,
                                  .paropts = NULL) {

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          fscale.response_spct,
          range = range,
          f = f,
          target = target,
          set.scaled = set.scaled,
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
                                set.scaled = target == 1,
                                ...,
                                .parallel = FALSE,
                                .paropts = NULL) {

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          fscale.filter_spct,
          range = range,
          f = f,
          target = target,
          set.scaled = set.scaled,
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
                                   set.scaled = target == 1,
                                   ...,
                                   .parallel = FALSE,
                                   .paropts = NULL) {

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          fscale.reflector_spct,
          range = range,
          f = f,
          target = target,
          set.scaled = set.scaled,
          qty.out = qty.out,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn fscale
#'
#' @export
#'
fscale.solute_mspct <- function(x,
                                range = NULL,
                                f = "mean",
                                target = 1,
                                set.scaled = target == 1,
                                ...,
                                .parallel = FALSE,
                                .paropts = NULL) {

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          fscale.solute_spct,
          range = range,
          f = f,
          target = target,
          set.scaled = set.scaled,
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
                             set.scaled = target == 1,
                             ...,
                             .parallel = FALSE,
                             .paropts = NULL) {

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          fscale.raw_spct,
          range = range,
          f = f,
          target = target,
          set.scaled = set.scaled,
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
                             set.scaled = target == 1,
                             ...,
                             .parallel = FALSE,
                             .paropts = NULL) {

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          fscale.cps_spct,
          range = range,
          f = f,
          target = target,
          set.scaled = set.scaled,
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
                                 set.scaled = target == 1,
                                 col.names,
                                 ...,
                                 .parallel = FALSE,
                                 .paropts = NULL) {

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          fscale, # members can be heterogeneous
          range = range,
          f = f,
          target = target,
          set.scaled = set.scaled,
          col.names = col.names,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}


# PRIVATE -----------------------------------------------------------------

#' fscale a spectrum
#'
#' This function returns a spectral object of the same class as the one
#' supplied as argument but with the spectral data re-scaled.
#'
#' @param spct generic_spct The spectrum to be normalized
#' @param range an R object on which range() returns a vector of length 2, with
#'   min and max wavelengths (nm)
#' @param col.names character The name of the variable to fscale
#' @param f function A summary function to be applied to \code{spct} or the
#'   name of such a function as a \code{character} string.
#' @param set.scaled logical Flag indicating if the data is to be marked
#'   as "scaled" or not.
#' @param na.rm logical Indicating whether NA values should be stripped before
#'   the computation proceeds
#' @param ... other arguments passed to f() by name.
#'
#' @return a new object of the same class as \code{spct}.
#'
#' @keywords internal
#'
fscale_spct <- function(spct,
                        range,
                        col.names,
                        f, target,
                        set.scaled,
                        na.rm = FALSE,
                        ...) {
  # Skip checks for intermediate results
  prev_state <- disable_check_spct()
  on.exit(set_check_spct(prev_state), add = TRUE)

  stopifnot("Missing columns" = all(col.names %in% colnames(spct)),
            "Multiple spectra in long form" = getMultipleWl(spct) == 1L)

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
  if (na.rm) {
    tmp.spct <- na.omit(tmp.spct)
  }
  multipliers <- numeric(length(col.names))
  i <- 0L
  for (col in col.names) {
    i <- i + 1L
    # rescaling needed
    if (!is.null(f)) {
      if (is.character(f)) {
        if (f %in% c("mean", "average")) {
          summary.value <- average_spct(tmp.spct[, c("w.length", col)], ...)
        } else if (f %in% c("total", "integral")) {
          summary.value <- integrate_spct(tmp.spct[, c("w.length", col)], ...)
        } else {
          f.name <- f
          f <- try(match.fun(f.name), silent = TRUE)
          if (is.function(f)) {
            summary.value <- f(tmp.spct[, c("w.length", col)], ...)
            f <- f.name
          } else {
            warning("A function named '", f.name, "' not found. Skipping!")
            summary.value <- NA_real_
            f <- NA
          }
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
  if (set.scaled) {
    spct <- setScaled(spct, scaled = list(multiplier = multipliers,
                                          f = f,
                                          range = range,
                                          target = target,
                                          cols = col.names))
  } else {
    spct <- setScaled(spct, scaled = FALSE)
  }
  check_spct(spct, force = FALSE)
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
#'
#' @family rescaling functions
#'
#' @examples
#'
#' scaled.spct <- fscale(sun.spct)
#' is_scaled(sun.spct)
#' is_scaled(scaled.spct)
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
#' @param .force.list logical If \code{TRUE} always silently return a
#'   list, with \code{FALSE} encoded field \code{multiplier = 1}.
#'
#' @return logical
#'
#' @note if x is not a \code{filter_spct} object, \code{NA} is returned
#'
#' @export
#'
#' @family rescaling functions
#'
#' @examples
#'
#' scaled.spct <- fscale(sun.spct)
#' getScaled(scaled.spct)
#'
getScaled <- function(x,
                      .force.list = FALSE) {

  if (is.generic_spct(x) || is.summary_generic_spct(x)) {

    scaled <- attr(x, "scaled", exact = TRUE)

    if (is.null(scaled) || (!is.list(scaled) && all(is.na(scaled))) ||
        (is.logical(scaled) && !scaled)) {
      # need to handle objects created with old versions
      if (.force.list) {
        scaled <- list(multiplier = 1,
                       f = NA,
                       range = c(NA_real_, NA_real_),
                       target = NA_real_,
                       cols = NA_character_)
      } else {
        scaled <- FALSE
      }
    } else if (is.list(scaled)) {
      if (!"target" %in% names(scaled)) {
        # cater for objects scaled before version 0.9.12
        scaled[["target"]] <- 1
      }
      if (!"range" %in% names(scaled)) {
        # cater for objects scaled before version 0.9.12
        scaled[["range"]] <- c(NA_real_, NA_real_)
      }
      if (!"cols" %in% names(scaled)) {
        # cater for objects scaled before version 0.10.10
        scaled[["cols"]] <- NA_character_
      }
    }
  } else {
    warning("Method 'getScaled()' not implemented for class: ",
            class(x)[1])
    if (.force.list) {
      scaled <- list(multiplier = NA_real_,
                     f = NA,
                     range = c(NA_real_, NA_real_),
                     target = NA_real_,
                     cols = NA_character_)
    } else {
      scaled <- NA
    }
  }
  scaled
}

#' @rdname getScaled
#'
#' @export
#'
getScaling <- function(x) {
  getScaled(x, .force.list = TRUE)
}

#' Set the "scaled" attribute
#'
#' Function to write the "scaled" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object.
#' @param scaled logical with \code{FALSE} meaning that values are expressed in
#'   absolute physical units and \code{TRUE} meaning that relative units are
#'   used. If \code{NULL} the attribute is not modified.
#' @param ... currently ignored.
#'
#' @note if x is not a \code{generic_spct} object, x is not modified.
#'
#' @export
#' @family rescaling functions
#'
setScaled <- function(x, ...) UseMethod("setScaled")

#' @describeIn setScaled Default for generic function
#'
#' @export
#'
#' @return a new object of the same class as \code{x}.
#'
setScaled.default <- function(x, ...) {
  warning("'setScaled()' is not defined for objects of class '", class(x)[1], "'.")
  return(x)
}

#' @describeIn setScaled Specialization for generic_spct
#'
#' @export
#'
#' @return a new object of the same class as \code{x}.
#'
setScaled.generic_spct <- function(x, ..., scaled = FALSE) {
  name <- substitute(x)
  if (!is.null(scaled)) {
    attr(x, "scaled") <- scaled
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' @describeIn setScaled Specialization for summary_generic_spct
#'
#' @export
#'
#' @return a new object of the same class as \code{x}.
#'
setScaled.summary_generic_spct <- setScaled.generic_spct

#' @describeIn setScaled Specialization for generic_mspct
#'
#' @export
#'
#' @return a new object of the same class as \code{x}.
#'
setScaled.generic_mspct <- function(x, ..., scaled = FALSE) {
  name <- substitute(x)
  z <- msmsply(x, setScaled.generic_spct, scaled = scaled)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, z, parent.frame(), inherits = TRUE)
  }
  invisible(z)
}
