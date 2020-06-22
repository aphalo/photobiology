
# normalize ---------------------------------------------------------------

#' Normalize spectral data
#'
#' This method returns a spectral object of the same class as the one
#' supplied as argument but with the spectral data normalized to 1.0 at a
#' specific wavelength.
#'
#' @details By default normalization is done based on the maximum of the
#'   spectral data. It is possible to also do the normalization based on a
#'   user-supplied wavelength expressed in nanometres or a wavelength retrieved
#'   using an arbitrary R function applied to the spectrum. By default the
#'   function is applied to the whole spectrum, but by passing a range of
#'   wavelengths as input, the search can be limited to a region of interest
#'   within the spectrum.
#'
#' @param x An R object
#' @param ... not used in current version
#'
#' @return A copy of \code{x}, with spectral data values normalized to one for
#' the criterion specified by the argument passed to \code{norm}.
#'
#' @examples
#'
#' normalize(sun.spct)
#' normalise(sun.spct) # equivalent
#'
#' normalize(sun.spct, norm = "max")
#' normalize(sun.spct, norm = 400)
#'
#' @export
#'
#' @family rescaling functions
#'
normalize <- function(x, ...) UseMethod("normalize")

#' @rdname normalize
#'
#' @note \code{normalise()} is a synonym for this \code{normalize()} method.
#'
#' @export
#'
normalise <- normalize

#' @describeIn normalize Default for generic function
#'
#' @export
#'
normalize.default <- function(x, ...) {
  warning("'normalize' is not defined for objects of class '", class(x)[1], "'.")
  x
}

#' @describeIn normalize Normalize a \code{source_spct} object.
#'
#' @param range An R object on which \code{range()} returns a numeric vector of
#'   length 2 with the limits of a range of wavelengths in nm, with min and max
#'   wavelengths (nm) used to set boundaries for search for normalization.
#' @param norm numeric Normalization wavelength (nm) or character string "max",
#'   or "min" for normalization at the corresponding wavelength.
#' @param unit.out character Allowed values "energy", and "photon",
#'   or its alias "quantum"
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before calculating the summary (e.g. "max") used for normalization.
#'
#' @note 1) By default if \code{x} contains one or more \code{NA} values and the
#'   normalization is based on a summary quantity, the returned spectrum will
#'   contain only \code{NA} values. If \code{na.rm == TRUE} then the summary
#'   quantity will be calculated after striping \code{NA} values, and only the
#'   values that were \code{NA} in \code{x} will be {NA} values in the returned
#'   spectrum.
#'
#' @return A copy of \code{x} with the values of the spectral quantity rescaled
#'   to 1 at the normalization wavelength. If the normalization wavelength is
#'   not already present in \code{x}, it is added by interpolation---i.e. the
#'   returned value may be one row longer than \code{x}.
#'
#' @export
#'
normalize.source_spct <- function(x,
                                  ...,
                                  range = NULL,
                                  norm = "max",
                                  unit.out = getOption("photobiology.radiation.unit",
                                                       default="energy"),
                                  na.rm = FALSE) {
  if (unit.out == "energy") {
    return(normalize_spct(spct = q2e(x, action = "replace"),
                          range = range,
                          norm = norm,
                          col.names = "s.e.irrad",
                          na.rm = na.rm))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(normalize_spct(spct = e2q(x, action = "replace"),
                          range = range,
                          norm = norm,
                          col.names = "s.q.irrad",
                          na.rm = na.rm))
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
                                  unit.out = getOption("photobiology.radiation.unit",
                                                       default="energy"),
                                  na.rm = FALSE) {
  if (unit.out == "energy") {
    return(normalize_spct(spct = q2e(x, action = "replace"),
                          range = range,
                          norm = norm,
                          col.names = "s.e.response",
                          na.rm = na.rm))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(normalize_spct(spct = e2q(x, action = "replace"),
                          range = range,
                          norm = norm,
                          col.names = "s.q.response",
                          na.rm = na.rm))
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
normalize.filter_spct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           qty.out = getOption("photobiology.filter.qty",
                               default = "transmittance"),
           na.rm = FALSE) {
    if (qty.out == "transmittance") {
      return(normalize_spct(spct = A2T(x, action = "replace"),
                            range = range,
                            norm = norm,
                            col.names = "Tfr",
                            na.rm = na.rm))
    } else if (qty.out == "absorbance") {
      return(normalize_spct(spct = T2A(x, action = "replace"),
                            range = range,
                            norm = norm,
                            col.names = "A",
                            na.rm = na.rm))
    } else if (qty.out == "absorptance") {
      return(normalize_spct(spct = T2Afr(x, action = "replace"),
                            range = range,
                            norm = norm,
                            col.names = "Afr",
                            na.rm = na.rm))
    } else {
      stop("'qty.out ", qty.out, " is unknown")
    }
  }

#' @describeIn normalize Normalize a reflector spectrum.
#'
#' @export
#'
normalize.reflector_spct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           qty.out = NULL,
           na.rm = FALSE) {
    normalize_spct(spct = x,
                   range = range,
                   norm = norm,
                   col.names = "Rfr",
                   na.rm = na.rm)
  }

#' @describeIn normalize Normalize a raw spectrum.
#'
#' @export
#'
normalize.raw_spct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           na.rm = FALSE) {
    normalize_spct(spct = x,
                   range = range,
                   norm = norm,
                   col.names = grep("^counts", names(x), value = TRUE),
                   na.rm = na.rm)
  }

#' @describeIn normalize Normalize a cps spectrum.
#'
#' @export
#'
normalize.cps_spct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           na.rm = FALSE) {
    normalize_spct(spct = x,
                   range = range,
                   norm = norm,
                   col.names = grep("^cps", names(x), value = TRUE),
                   na.rm = na.rm)
  }

#' @describeIn normalize Normalize a raw spectrum.
#'
#' @param col.names character vector containing the names of columns or
#'   variables to which to apply the normalization.
#'
#' @export
#'
normalize.generic_spct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           col.names,
           na.rm = FALSE) {

    normalize_spct(spct = x,
                   range = range,
                   norm = norm,
                   col.names = col.names,
                   na.rm = na.rm)
  }

# collections of spectra --------------------------------------------------


#' @describeIn normalize Normalize the members of a source_mspct object.
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
normalize.source_mspct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           na.rm = FALSE,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            normalize,
            range = range,
            norm = norm,
            unit.out = unit.out,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)

  }

#' @describeIn normalize Normalize the members of a response_mspct object.
#'
#' @export
#'
normalize.response_mspct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           na.rm = FALSE,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            normalize,
            range = range,
            norm = norm,
            unit.out = unit.out,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)

  }

#' @describeIn normalize Normalize the members of a filter_mspct object.
#'
#' @export
#'
normalize.filter_mspct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           qty.out = getOption("photobiology.filter.qty",
                               default = "transmittance"),
           na.rm = FALSE,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            normalize,
            range = range,
            norm = norm,
            qty.out = qty.out,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)

  }

#' @describeIn normalize Normalize the members of a reflector_mspct object.
#'
#' @export
#'
normalize.reflector_mspct <- function(x,
                                      ...,
                                      range = x,
                                      norm = "max",
                                      qty.out = NULL,
                                      na.rm = FALSE,
                                      .parallel = FALSE,
                                      .paropts = NULL) {
  msmsply(x,
          normalize,
          range = range,
          norm = norm,
          qty.out = qty.out,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)

}

#' @describeIn normalize Normalize the members of a raw_mspct object.
#'
#' @export
#'
normalize.raw_mspct <- function(x,
                                ...,
                                range = x,
                                norm = "max",
                                na.rm = FALSE,
                                .parallel = FALSE,
                                .paropts = NULL) {
  msmsply(x,
          normalize,
          range = range,
          norm = norm,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn normalize Normalize the members of a cps_mspct object.
#'
#' @export
#'
normalize.cps_mspct <- function(x,
                                ...,
                                range = x,
                                norm = "max",
                                na.rm = FALSE,
                                .parallel = FALSE,
                                .paropts = NULL) {
  msmsply(x,
          normalize,
          range = range,
          norm = norm,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

# PRIVATE -----------------------------------------------------------------

#' @keywords internal
#'
normalize_spct <- function(spct, range, norm, col.names, na.rm) {
  stopifnot(is.generic_spct(spct), !is.null(col.names),
            col.names %in% names(spct))
  if (getMultipleWl(spct) != 1L) {
    warning("Object contains data for ",
            getMultipleWl(spct), " spectra; skipping normalization")
    return(spct)
  }
  # normalization will wipe out any existing scaling
  if (is_scaled(spct)) {
    setScaled(spct, scaled = FALSE)
  }

  if (is_normalized(spct)) {
    setNormalized(spct, norm = FALSE)
  }

  if (na.rm) {
    x <- na.omit(spct)
  } else {
    x <- spct
  }

  if (is.null(range) || all(is.na(range))) {
    range <- range(spct)
  } else {
    x <- trim_wl(x, range)
  }
  stopifnot(nrow(x) > 2) # too short a slice

  # rescaling needed
  if (!is.null(norm)) {
    for (col in col.names) {
      if (is.character(norm)) {
        if (norm %in% c("max", "maximum")) {
          idx <- which.max(x[[col]])
        } else if (norm %in% c("min", "minimum")) {
          idx <- which.min(x[[col]])
        } else {
          warning("Invalid character '", norm, "'value in 'norm'")
          idx <- NA
        }
        scale.factor <- 1 / x[idx, col, drop = TRUE]
        norm <- x[idx, "w.length", drop = TRUE]
      } else if (is.numeric(norm)) {
        if (norm >= range[1] && norm <= range[2]) {
          tmp.spct <- spct[ , c("w.length", col)]
          class(tmp.spct) <- class(spct)
          scale.factor <- 1 /
            interpolate_spct(spct = tmp.spct, w.length.out = norm)[ , eval(col)]
        } else {
          warning("'norm = ", norm, "' value outside spectral data range of ",
                  round(min(tmp.spct), 1), " to ", round(max(tmp.spct), 1), " (nm)")
          scale.factor <- NA
        }
      } else {
        stop("'norm' should be numeric or character")
      }
      spct[[col]] <- spct[ , col, drop = TRUE] * scale.factor
    }
  } else {
    return(spct)
  }
  spct <- setNormalized(spct, norm)
  spct
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
  if (!is.generic_spct(x) && !is.summary_generic_spct(x)) {
    return(NA)
  }
  spct.attr <- attr(x, "normalized", exact = TRUE)
  !is.null(spct.attr) && as.logical(spct.attr)
}

#' @rdname is_normalized
#'
#' @note \code{is_normalised()} is a synonym for this \code{is_normalized()}
#'   method.
#'
#' @export
#'
is_normalised <- is_normalized

# getNormalized -----------------------------------------------------------

#' Get the "normalized" attribute
#'
#' Function to read the "normalized" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#' @param .force.numeric logical If \code{TRUE} always silently return a
#'   numeric value, with \code{FALSE} encoded as zero, and character values
#'   as \code{NA}.
#'
#' @return numeric or logical (possibly character for objects created with
#'   earlier versions).
#'
#' @note if x is not a \code{generic_spct} object, \code{NA} is returned
#'
#' @export
#' @family rescaling functions
#'
getNormalized <- function(x,
                          .force.numeric = FALSE) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    normalized <- attr(x, "normalized", exact = TRUE)
    if (is.null(normalized) || is.na(normalized)) {
      # need to handle objects created with old versions
      normalized <- FALSE
    }
  } else {
    warning("Method 'getNormalized()' not implemented for class: ",
            class(x)[1])
    normalized <- NA
  }
  if (.force.numeric) {
    suppressWarnings(as.numeric(normalized[[1]]))
  } else {
    normalized[[1]]
  }
}

#' @rdname getNormalized
#'
#' @note \code{getNormalised()} is a synonym for this \code{getNormalized()}
#'   method.
#'
#' @export
#'
getNormalised <- getNormalized

#' Set the "normalized" attribute
#'
#' Function to write the "normalized" attribute of an existing generic_spct
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
  if ((is.generic_spct(x) || is.summary_generic_spct(x)) &&
      (is.na(norm) || is.numeric(norm) || is.logical(norm))) {
    attr(x, "normalized") <- norm
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' @rdname setNormalized
#'
#' @note \code{setNormalised()} is a synonym for this \code{setNormalized()}
#'   method.
#'
#' @export
#'
setNormalised <- setNormalized
