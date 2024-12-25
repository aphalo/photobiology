
# normalize ---------------------------------------------------------------

#' Normalize spectral data
#'
#' This method returns a spectral object of the same class as the one supplied
#' as argument but with the spectral data normalized to 1.0 at a specific
#' wavelength.  When the object contains multiple spectra, the normalisation is
#' applied to each spectrum individually.
#'
#' @details By default normalization is done based on the maximum of the
#'   spectral data. It is possible to also do the normalization based on a
#'   user-supplied wavelength expressed in nanometres or the minimum. An
#'   existing normalization can be updated for a different unit of expression or
#'   after a conversion to a related spectral quantity.
#'
#'   By default the function is applied to the whole spectrum, but by passing a
#'   range of wavelengths as input, the search, e.g., for the maximum, can be
#'   limited to a range of wavelengths of interest instead of the whole
#'   spectrum.
#'
#'   In 'photobiology' (>= 0.10.8) detailed information about the normalization
#'   is stored in an attribute. In 'photobiology' (>= 0.10.10)
#'   applying a new normalization to an already normalized spectrum recomputes
#'   the multiplier factors stored in the attributes whenever possible. This
#'   ensures that the returned object is identical, except for possible
#'   accumulated loss of precision due to floating-point arithmetic,
#'   independently of the previous application of a different normalization.
#'
#' @note When the spectrum passed as argument to \code{x} had been previously
#'   scaled, in 'photobiology' (<= 0.10.9) the scaling attribute was always
#'   removed and no normalization factors returned. In 'photobiology'
#'   (>= 0.10.10) scaling information can be preserved by passing
#'   \code{keep.scaling = TRUE}.
#'
#'   By default if \code{x} contains one or more \code{NA} values and the
#'   normalization is based on a summary quantity, the returned spectrum will
#'   contain only \code{NA} values. If \code{na.rm == TRUE} then the summary
#'   quantity will be calculated after striping \code{NA} values, and only the
#'   values that were \code{NA} in \code{x} will be \code{NA} values in the
#'   returned spectrum.
#'
#'   When a numeric value is passed as argument to keep.scaling, the scaling
#'   uses \code{f = "total"} or \code{f = "mean"} depending on the class of
#'   \code{x}. Prescaling is only occasionally needed.
#'
#'   Method \code{normalize} is implemented for \code{solute_spct} objects but
#'   as the spectral data stored in them are a description of an intensive
#'   property of a substance, normalization is unlikely to useful. To represent
#'   solutions of specific concentrations of solutes, \code{filter_spct} objects
#'   should be used instead.
#'
#' @param x An R object
#' @param ... not used in current version
#'
#' @return A copy of the object passed as argument to \code{x} with the values
#'   of the spectral quantity rescaled to 1 at the normalization wavelength. If
#'   the normalization wavelength is not already present in \code{x}, it is
#'   added by interpolation---i.e. the returned value may be one row longer than
#'   \code{x}. Attributes \code{normalized} and \code{normalization} are set to
#'   keep a log of the computations applied.
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
#'   or "min" for normalization at the corresponding wavelength, "update" to
#'   update the normalization after modifying units of expression, quantity
#'   or range but respecting the previously used criterion, "undo" to revert
#'   an existing normalization or "skip" to force return of \code{x} unchanged.
#' @param unit.out No longer supported and is ignored with a warning.
#' @param keep.scaling logical or numeric Flag to indicate if any existing
#'   scaling should be preserved or not. The default, \code{FALSE}, preserves
#'   the behaviour of versions (<= 0.10.9). If numeric, the spectrum is scaled
#'   to this value before normalization and marked as not scaled.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before calculating the summary (e.g. "max") used for normalization.
#'
#' @export
#'
normalize.source_spct <- function(x,
                                  ...,
                                  range = NULL,
                                  norm = "max",
                                  unit.out = NA,
                                  keep.scaling = FALSE,
                                  na.rm = FALSE) {
  if (!is.na(unit.out)) {
    warning("Argument 'unit.out' is no longer supported and is ignored.")
  }

  if (getMultipleWl(x) > 1L) {
    # brute force and slow approach, unsuitable for long time series
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      normalize.source_mspct(x = mspct,
                             range = range,
                             norm = norm,
                             keep.scaling = keep.scaling,
                             na.rm = na.rm,
                             ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  col.names <- intersect(colnames(x), c("s.e.irrad", "s.q.irrad"))
  return(normalize_spct(spct = x,
                        range = range,
                        norm = norm,
                        col.names = col.names,
                        keep.scaling = keep.scaling,
                        na.rm = na.rm,
                        ...))
}

#' @describeIn normalize Normalize a response spectrum.
#'
#' @export
#'
normalize.response_spct <- function(x,
                                    ...,
                                    range = NULL,
                                    norm = "max",
                                    unit.out = NA,
                                    keep.scaling = FALSE,
                                    na.rm = FALSE) {
  if (!is.na(unit.out)) {
    warning("Argument 'unit.out' is no longer supported and is ignored.")
  }

  if (getMultipleWl(x) > 1L) {
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      normalize.response_mspct(x = mspct,
                               range = range,
                               norm = norm,
                               keep.scaling = keep.scaling,
                               na.rm = na.rm,
                               ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  col.names <- intersect(colnames(x), c("s.e.response", "s.q.response"))
  return(normalize_spct(spct = x,
                        range = range,
                        norm = norm,
                        col.names = col.names,
                        keep.scaling = keep.scaling,
                        na.rm = na.rm,
                        ...))
}

#' @describeIn normalize Normalize a filter spectrum.
#'
#' @param qty.out No longer supported and is ignored with a warning..
#'
#' @export
#'
normalize.filter_spct <- function(x,
                                  ...,
                                  range = NULL,
                                  norm = "max",
                                  qty.out = NA,
                                  keep.scaling = FALSE,
                                  na.rm = FALSE) {
  if (!is.na(qty.out)) {
    warning("Argument 'qty.out' is no longer supported and is ignored.")
    }

    if (getMultipleWl(x) > 1L) {
      mspct <- subset2mspct(x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      mspct <-
        normalize.filter_mspct(x = mspct,
                               range = range,
                               norm = norm,
                               keep.scaling = keep.scaling,
                               na.rm = na.rm,
                               ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    col.names <- intersect(colnames(x), c("Tfr", "A", "Afr"))
    return(normalize_spct(spct = x,
                          range = range,
                          norm = norm,
                          col.names = col.names,
                          keep.scaling = keep.scaling,
                          na.rm = na.rm,
                          ...))
  }

#' @describeIn normalize Normalize a reflector spectrum.
#'
#' @export
#'
normalize.reflector_spct <- function(x,
                                     ...,
                                     range = NULL,
                                     norm = "max",
                                     qty.out = NA,
                                     keep.scaling = FALSE,
                                     na.rm = FALSE) {
    if (!is.na(qty.out)) {
      warning("Argument 'qty.out' is no longer supported and is ignored.")
    }

    if (getMultipleWl(x) > 1L) {
      mspct <- subset2mspct(x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      mspct <-
        normalize.reflector_mspct(x = mspct,
                                  range = range,
                                  norm = norm,
                                  keep.scaling = keep.scaling,
                                  na.rm = na.rm,
                                  ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    normalize_spct(spct = x,
                   range = range,
                   norm = norm,
                   col.names = "Rfr",
                   keep.scaling = keep.scaling,
                   na.rm = na.rm,
                   ...)
  }

#' @describeIn normalize Normalize a solute spectrum.
#'
#' @export
#'
normalize.solute_spct <- function(x,
                                  ...,
                                  range = NULL,
                                  norm = "max",
                                  qty.out = NA,
                                  keep.scaling = FALSE,
                                  na.rm = FALSE) {
  if (!is.na(qty.out)) {
    warning("Argument 'qty.out' is no longer supported and is ignored.")
  }

  # for consistency use qty.out parameter and add support!!!
  if (getMultipleWl(x) > 1L) {
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      normalize.solute_mspct(x = mspct,
                             range = range,
                             norm = norm,
                             keep.scaling = keep.scaling,
                             na.rm = na.rm,
                             ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  col.names <- intersect(c("K.mole", "K.mass"), colnames(x))
  normalize_spct(spct = x,
                 range = range,
                 norm = norm,
                 col.names = col.names,
                 keep.scaling = keep.scaling,
                 na.rm = na.rm,
                 ...)
}

#' @describeIn normalize Normalize a raw spectrum.
#'
#' @export
#'
normalize.raw_spct <- function(x,
                               ...,
                               range = NULL,
                               norm = "max",
                               keep.scaling = FALSE,
                               na.rm = FALSE) {
  if (getMultipleWl(x) > 1L) {
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      normalize.raw_mspct(x = mspct,
                          range = range,
                          norm = norm,
                          keep.scaling = keep.scaling,
                          na.rm = na.rm,
                          ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  col.names <- grep("^counts", names(x), value = TRUE)
  normalize_spct(spct = x,
                 range = range,
                 norm = norm,
                 col.names = col.names,
                 keep.scaling = keep.scaling,
                 na.rm = na.rm,
                 ...)
}

#' @describeIn normalize Normalize a cps spectrum.
#'
#' @export
#'
normalize.cps_spct <- function(x,
                               ...,
                               range = NULL,
                               norm = "max",
                               keep.scaling = FALSE,
                               na.rm = FALSE) {
  if (getMultipleWl(x) > 1L) {
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      normalize.cps_mspct(x = mspct,
                          range = range,
                          norm = norm,
                          keep.scaling = keep.scaling,
                          na.rm = na.rm,
                          ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  col.names <- grep("^cps", names(x), value = TRUE)
  normalize_spct(spct = x,
                 range = range,
                 norm = norm,
                 col.names = grep("^cps", names(x), value = TRUE),
                 keep.scaling = keep.scaling,
                 na.rm = na.rm,
                 ...)
}

#' @describeIn normalize Normalize a raw spectrum.
#'
#' @param col.names character vector containing the names of columns or
#'   variables. Columns in \code{x} matching the names in \code{col.names} are
#'   normalized, other columns are returned unchanged.
#'
#' @export
#'
normalize.generic_spct <- function(x,
                                   ...,
                                   range = NULL,
                                   norm = "max",
                                   col.names,
                                   keep.scaling = FALSE,
                                   na.rm = FALSE) {
  if (getMultipleWl(x) > 1L) {
    mspct <- subset2mspct(x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    mspct <-
      normalize.generic_mspct(x = mspct,
                              range = range,
                              norm = norm,
                              col.names = col.names,
                              keep.scaling = keep.scaling,
                              na.rm = na.rm,
                              ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }
  col.names <- intersect(col.names, colnames(x))
  if (length(col.names)) {
    normalize_spct(spct = x,
                   range = range,
                   norm = norm,
                   col.names = col.names,
                   keep.scaling = keep.scaling,
                   na.rm = na.rm,
                   ...)
  } else {
    message("No columns to normalize.")
    x
  }
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
normalize.source_mspct <- function(x,
                                   ...,
                                   range = NULL,
                                   norm = "max",
                                   unit.out = NA,
                                   keep.scaling = FALSE,
                                   na.rm = FALSE,
                                   .parallel = FALSE,
                                   .paropts = NULL) {
  if (!is.na(unit.out)) {
    warning("Argument 'unit.out' is no longer supported and is ignored.")
  }

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          normalize.source_spct,
          range = range,
          norm = norm,
          keep.scaling = keep.scaling,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)

}

#' @describeIn normalize Normalize the members of a response_mspct object.
#'
#' @export
#'
normalize.response_mspct <- function(x,
                                     ...,
                                     range = NULL,
                                     norm = "max",
                                     unit.out = NA,
                                     keep.scaling = FALSE,
                                     na.rm = FALSE,
                                     .parallel = FALSE,
                                     .paropts = NULL) {
  if (!is.na(unit.out)) {
    warning("Argument 'unit.out' is no longer supported and is ignored.")
  }

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          normalize.response_spct,
          range = range,
          norm = norm,
          keep.scaling = keep.scaling,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)

}

#' @describeIn normalize Normalize the members of a filter_mspct object.
#'
#' @export
#'
normalize.filter_mspct <- function(x,
                                   ...,
                                   range = NULL,
                                   norm = "max",
                                   qty.out = NA,
                                   keep.scaling = FALSE,
                                   na.rm = FALSE,
                                   .parallel = FALSE,
                                   .paropts = NULL) {

  if (!is.na(qty.out)) {
    warning("Argument 'qty.out' is no longer supported and is ignored.")
  }

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          normalize.filter_spct,
          range = range,
          norm = norm,
          keep.scaling = keep.scaling,
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
                                      qty.out = NA,
                                      keep.scaling = FALSE,
                                      na.rm = FALSE,
                                      .parallel = FALSE,
                                      .paropts = NULL) {
  if (!is.na(qty.out)) {
    warning("Argument 'qty.out' is no longer supported and is ignored.")
  }

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          normalize.reflector_spct,
          range = range,
          norm = norm,
          keep.scaling = keep.scaling,
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
                                keep.scaling = FALSE,
                                na.rm = FALSE,
                                .parallel = FALSE,
                                .paropts = NULL) {

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          normalize.raw_spct,
          range = range,
          norm = norm,
          keep.scaling = keep.scaling,
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
                                keep.scaling = FALSE,
                                na.rm = FALSE,
                                .parallel = FALSE,
                                .paropts = NULL) {

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          normalize.cps_spct,
          range = range,
          norm = norm,
          keep.scaling = keep.scaling,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn normalize Normalize the members of a solute_mspct object.
#'
#' @export
#'
normalize.solute_mspct <- function(x,
                                   ...,
                                   range = x,
                                   norm = "max",
                                   qty.out = NA,
                                   keep.scaling = FALSE,
                                   na.rm = FALSE,
                                   .parallel = FALSE,
                                   .paropts = NULL) {
  if (!is.na(qty.out)) {
    warning("Argument 'qty.out' is no longer supported and is ignored.")
  }

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          normalize.solute_spct,
          range = range,
          norm = norm,
          keep.scaling = keep.scaling,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn normalize Normalize the members of a solute_mspct object.
#'
#' @export
#'
normalize.generic_mspct <- function(x,
                                    ...,
                                    range = NULL,
                                    norm = "max",
                                    col.names,
                                    keep.scaling = FALSE,
                                    na.rm = FALSE,
                                    .parallel = FALSE,
                                    .paropts = NULL) {

  if (!length(x)) return(x) # class of x in no case changes

  msmsply(x,
          normalize, # members can be heterogeneous
          range = range,
          norm = norm,
          col.names = col.names,
          keep.scaling = keep.scaling,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

# PRIVATE -----------------------------------------------------------------

#' @keywords internal
#'
normalize_spct <- function(spct,
                           range,
                           norm,
                           col.names,
                           na.rm,
                           keep.scaling,
                           ...) {
  stopifnot(is.generic_spct(spct))

  # if 'norm' is a character vector, we use the first element
  # thus, all columns always get the same type of normalization
  if (is.character(norm) && length(norm) > 1) {
    if (length(unique(norm)) > 1) {
      warning("Multiple 'norm' values supplied by name. Using the first one: ",
              norm[1], ".")
    }
    morm <- norm[1]
  }

  # handle "skip" early so that long-form multiple spectra or missing columns
  # do not trigger errors
  if (!length(norm) ||
      any(is.na(norm)) ||
      norm[1] == "skip" ||
      (norm[1] == "update" && !is_normalized(spct))) {
    return(spct)
  } else {
    norm <- rep_len(norm, length(col.names))
  }

  stopifnot("Missing columns" = all(col.names %in% colnames(spct)),
            "Multiple spectra in long form" = getMultipleWl(spct) == 1L)

  updating <- is_normalized(spct)

  if (updating) {
    # we retrieve the existing normalization data
    old.normalization.ls <- getNormalization(spct)

    required.fields <- c("norm.type", "norm.wl", "norm.cols", "norm.range")
    has.normalization.metadata <-
      !any(is.na(unlist(old.normalization.ls[required.fields])))

    if (!has.normalization.metadata) {
      warning("Normalization not updated: action not supported for objects created with 'photobiology' (<= 0.10.9).")
      return(spct)
    } else {
      if (norm[1] != "undo") {
        norm <- old.normalization.ls$norm.type
        if (norm[1] == "wavelength") {
          norm <- old.normalization.ls$norm.wl
        }
        range <- old.normalization.ls$norm.range
      }
      # remove the old normalization
      spct <- denormalize_spct(spct, wipe.away = FALSE)
      if (norm[1] == "undo") {
        return(spct)
      }
    }
  } else if (norm[1] == "update") {
    # not normalized, nothing to update
    return(spct)
  }

  # we do not remove NA's from the returned object, only for computation
  if (na.rm) {
    x <- na.omit(spct)
  } else {
    x <- spct
  }
  # set or check 'range'
  if (is.null(range) || all(is.na(range))) {
    range <- wl_range(x, na.rm = TRUE)
  } else {
    x <- trim_wl(x, range)
    range <- wl_range(x, na.rm = TRUE) # if range was broader x is not expanded
  }
  stopifnot("Wavelength 'range' is too narrow" = nrow(x) > 2) # too short a slice

  norm.arg <- norm
  # normalization will wipe out any existing scaling except for its effect
  # on the computed factors.
  if (is_scaled(x) && !keep.scaling) {
    # Only behaviour in <= 0.10.9
    # remove scaling metadata and do not save norm.factors
    scale.is.dirty <- TRUE
    setScaled(spct, scaled = FALSE)
  } else {
    # retain scaling metadata and save norm.factors
    scale.is.dirty <- FALSE
  }

  # normalization of one or more columns
  scale.factors <- numeric(0)
  norm.wls <- numeric(0)
  for (i in seq_along(col.names)) {
    col <- col.names[i]
    if (is.character(norm[i])) {
      if (norm[i] %in% c("max", "maximum")) {
        idx <- which.max(x[[col]])
      } else if (norm[i] %in% c("min", "minimum")) {
        idx <- which.min(x[[col]])
      } else {
        warning("Invalid 'norm' value: '", norm[i], "'")
        idx <- NA
      }
      scale.factor <- 1 / x[idx, col, drop = TRUE]
      norm.wl <- x[idx, "w.length", drop = TRUE]
    } else if (is.numeric(norm)) {
      if (norm[i] >= range[1] && norm[i] <= range[2]) {
        norm.wl <- norm[i]
        tmp.spct <- spct[ , c("w.length", col)]
        class(tmp.spct) <- class(spct)
        scale.factor <- 1 /
          interpolate_spct(spct = tmp.spct, w.length.out = norm.wl)[ , eval(col)]
      } else {
        warning("'norm = ", norm[i], "' value(s) outside spectral data range of ",
                round(min(tmp.spct), 1), " to ", round(max(tmp.spct), 1), " (nm)")
        scale.factor <- NA
      }
    } else {
      stop("'norm' should be numeric or character")
    }
    spct[[col]] <- spct[ , col, drop = TRUE] * scale.factor
    scale.factors <- c(scale.factors, scale.factor)
    norm.wls <- c(norm.wls, norm.wl)
  }

  z <- setNormalized(spct,
                     norm = norm.wls,
                     norm.type =
                       if (is.character(norm.arg)) {
                         norm.arg
                       } else if (is.numeric(norm.arg)) {
                         "wavelength"
                       },
                     norm.factors =
                       if (scale.is.dirty) {
                         rep(NA_real_, length(col.names))
                       } else {
                         scale.factors
                       },
                     norm.cols = col.names,
                     norm.range = range)
  z # setNormalized makes its returned value invisible
}

#' Remove a previously applied normalization
#'
#' @param spct An object of one of the classes derived from \code{geberic_spct}
#'   containing data for a single spectrum (\code{getMultipleWl(spct) == 1}).
#' @param wipe.away logical If \code{TRUE} the normalization metadata is removed
#'   without undoing the effect of the normalization.
#'
#' @return A modified copy of \code{spct} if it was previously normalized or
#' \code{spct} unchanged, otherwise.
#'
#' @keywords internal
#'
denormalize_spct <- function(spct, wipe.away = FALSE) {
  if (!is_normalized(spct)) {
    return(spct)
  }

  if (wipe.away) {
    message("Removing normalization metadata keeping normalization!")
    attr(spct, "normalized") <- FALSE
    attr(spct, "normalization") <- NULL
    return(spct)
  }

  old.normalization.ls <- getNormalization(spct)
  required.fields <-
    c("norm.factors", "norm.cols")
  has.normalization.metadata <-
    !any(is.na(unlist(old.normalization.ls[required.fields])))
  norm.ls.idx <- which(old.normalization.ls$norm.cols %in% colnames(spct))

  if (has.normalization.metadata) {
    stopifnot("Missing columns" = length(norm.ls.idx) > 0L,
              "Multiple spectra in long form" = getMultipleWl(spct) == 1L)
    # earlier bug lead to not saving all norm.factors for multiple columns
    if (length(old.normalization.ls$norm.cols) !=
        length(old.normalization.ls$norm.factors)) {
      warning("Normalization metadata incomplete, denormalization not possible.")
      return(spct)
    }
    for (i in seq_along(old.normalization.ls$norm.cols)) {
      if (!i %in% norm.ls.idx) {
        next()
      }
      col.name.i <- old.normalization.ls$norm.cols[i]
      norm.factor.i <- old.normalization.ls$norm.factors[i]
      spct[[col.name.i]] <- spct[[col.name.i]] / norm.factor.i
    }
    attr(spct, "normalized") <- FALSE
    attr(spct, "normalization") <- NULL
    return(spct)
  } else {
    stop("Normalization metadata missing, denormalization not possible.")
  }
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
#' @return A \code{logical} value indicating if \code{x} is normalized or not,
#'   for collections of spectra, a named list with \code{logicals} as members.
#'   If \code{x} is not a \code{generic_spct} or \code{generic_mspct} object the
#'   value returned is \code{NA}.
#'
#' @export
#' @family rescaling functions
#'
is_normalized <- function(x) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    spct.attr <- attr(x, "normalized", exact = TRUE)
    # in some versions a logical was used, but later the normalization wavelength
    # in old versions the attribute was set only when normalization was applied
    stopifnot(is.null(spct.attr) || is.numeric(spct.attr) || is.logical(spct.attr))
    return(!is.null(spct.attr) && as.logical(spct.attr))
  } else if (is.generic_mspct(x)) {
    return(mslply(x, is_normalized))
  } else if (is.waveband(x)) {
    return(!is.na(normalization(x)))
  } else {
    warning("Method 'is_normalized()' not implemented for class: '", class(x)[1], "'.")
    return(NA)
  }
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

#' Query the "normalized" and "normalization" attributes
#'
#' Functions to read the "normalized" and "normalization" attributes of an
#' existing generic_spct object.
#'
#' @param x a generic_spct object.
#' @param .force.numeric logical If \code{TRUE} always silently return a
#'   numeric value, with \code{FALSE} encoded as zero, and character values
#'   as \code{NA}.
#'
#' @details Spectral data that has been normalized needs to be used diffferently
#'   in computations than data expresed in original units. These two functions
#'   make it possible to query if data stored in an object of class
#'   \code{generic_spct} or of a derived class contains data expressed in
#'   physical units or normalized. In the later case, it is possible to also
#'   query how the normalization was done.
#'
#' @return \code{getNormalized()} returns numeric or logical (possibly character
#'   for objects created with earlier versions); for collections of spectra, a
#'   named list, with one member for each spectrum. If \code{x} is not a
#'   \code{generic_spct} object, \code{NA} or a list with fields set to NAs is
#'   returned. Objects created with versions of package 'photobiology' earlier
#'   than 0.10.8 are lacking the detailed normalization metadata.
#'
#' @export
#'
#' @examples
#'
#' getNormalized(sun.spct)
#' getNormalization(sun.spct)
#'
#' sun_norm.spct <- normalize(sun.spct)
#'
#' getNormalized(sun_norm.spct)
#' getNormalization(sun_norm.spct)
#'
#' getNormalization(e2q(sun_norm.spct))
#'
#' gel_norm.spct <- normalize(yellow_gel.spct)
#'
#' getNormalized(gel_norm.spct)
#' getNormalization(gel_norm.spct)
#'
#' # getNormalization(T2Afr(gel_norm.spct))
#' getNormalization(any2A(gel_norm.spct))
#'
#' @family rescaling functions
#'
getNormalized <- function(x,
                          .force.numeric = FALSE) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    normalized <- attr(x, "normalized", exact = TRUE)
    if (is.null(normalized) || is.na(normalized)) {
      # need to handle objects created with very old versions
      normalized <- FALSE
    }
  } else if (is.generic_mspct(x)) {
    return(mslply(x, getNormalized, .force.numeric = .force.numeric))
  } else {
    warning("Method 'getNormalized()' not implemented for class: ",
            class(x)[1])
    normalized <- NA
  }
  if (.force.numeric) {
    suppressWarnings(as.numeric(normalized))
  } else {
    normalized
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

#' @rdname getNormalized
#'
#' @return \code{getNormalization()} returns a list with five fields: norm.type,
#'   norm.wl, norm.factors, norm.cols, norm.range. For collections of spectra, a
#'   named list of lists, with one member list for each member of the collection
#'   of spectra. See \code{\link{setNormalized}()} for the values stored in the
#'   fields.
#'
#' @export
#'
getNormalization <- function(x) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    if (is_normalized(x)) {
      # attribute in use >= 0.10.8
      normalization.list <- attr(x, "normalization", exact = TRUE)
      if (is.list(normalization.list)) {
        if (!exists("norm.range", normalization.list)) {
          # norm.range is missing 0.10.8 and 0.10.9
          normalization.list[["norm.range"]] <- rep(NA_real_, 2)
        }
        return(normalization.list)
      } else if (is.numeric(getNormalized(x, .force.numeric = FALSE))) {
        return(list(norm.type = NA_character_,
                    norm.wl = getNormalized(x),
                    norm.factors = NA_real_,
                    norm.cols = NA_character_,
                    norm.range = rep(NA_real_, 2))
        )
      }
    }
  } else if (is.generic_mspct(x)) {
    return(mslply(x, getNormalization))
  } else {
    warning("Method 'getNormalization()' not implemented for class: ",
            class(x)[1])
  }
  list(norm.type = NA_character_,
       norm.wl = NA_real_,
       norm.factors = NA_real_,
       norm.cols = NA_character_,
       norm.range = rep(NA_real_, 2))
}

#' @rdname getNormalized
#' @export
#'
getNormalisation <- getNormalization

#' Set the "normalized" and "normalization" attributes
#'
#' Function to write the "normalized" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object.
#' @param norm numeric (or logical) Normalization wavelength (nanometres).
#' @param norm.type character Type of normalization applied.
#' @param norm.factors numeric The scaling factor(s) so that dividing the spectral
#'   values by this factor reverts the normalization.
#' @param norm.cols character The name(s) of the data columns normalized.
#' @param norm.range numeric The wavelength range used for normalization (nm).
#' @param verbose logical Flag enabling or silencing informative warnings.
#'
#' @details This function \strong{is used internally}, although occasionally
#'   users may want to use it to "pretend" that spectral data have not been
#'   normalized. Use \code{\link{normalize}()} methods to apply a normalization
#'   and set the attributes accordingly. Function \code{setNormalized()} only
#'   sets the attributes that store the metadata corresponding to an already
#'   applied normalization. Thus a trace of the transformations applied to
#'   spectral data is kept, which currently is used to renormalize the spectra
#'   when the quantity used for expression is changed with a conversion
#'   function. It is also used in other packages like 'ggspectra' when
#'   generating automatically axis labels. If \code{x} is not a
#'   \code{generic_spct} object, \code{x} is not modified.
#'
#' @note Passing a \code{logical} as argument to \code{norm} is deprecated
#'   but accepted silently for backwards compatibility.
#'
#' @export
#' @family rescaling functions
#'
setNormalized <- function(x,
                          norm = FALSE,
                          norm.type = NA_character_,
                          norm.factors = NA_real_,
                          norm.cols = NA_character_,
                          norm.range = rep(NA_real_, 2),
                          verbose = getOption("verbose_as_default", default = FALSE)) {
  name <- substitute(x)

  stopifnot("'norm' must be numeric or logical, but it is not" =
              is.numeric(norm) || is.logical(norm))

  if (!length(norm.cols)) {
    norm.cols = NA_character_
  } else if (!all(is.na(norm.cols))) {
    norm.cols <- na.omit(norm.cols)
    norm.cols <- intersect(norm.cols, setdiff(colnames(x), "w.length"))
    stopifnot("'norm.cols' argument must be 'NA' or name(s) of columns" =
                length(norm.cols) > 0L)
  }

  if (is.logical(norm) && all(!norm)) {
    attr(x, "normalized") <- FALSE
    attr(x, "normalization") <- NULL
  } else if ((is.generic_spct(x) || is.summary_generic_spct(x)) &&
             (all(is.na(norm)) || all(is.numeric(norm)) || all(is.logical(norm)))) {
    attr(x, "normalized") <-
      if (length(norm) == 1L) {
        norm
      } else {
        TRUE
      }
    normalization.ls <- list(norm.type = norm.type,
                             norm.wl = if(is.numeric(norm)) norm else NA_real_,
                             norm.factors = norm.factors,
                             norm.cols = norm.cols,
                             norm.range = norm.range)
    if (verbose && norm && anyNA(normalization.ls, recursive = TRUE)) {
      message("\"normalized\" attribute set to TRUE, with missing ",
              paste(names(normalization.ls)[which(is.na(normalization.ls))], collapse = ". "),
              "data.")
    }
    attr(x, "normalization") <- normalization.ls
  } else {
    warning("Method 'setNormalization()' not implemented for class: ",
            class(x)[1])
    invisible(x) # return the object unchanged
  }
  # set by reference
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
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
