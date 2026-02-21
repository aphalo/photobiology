
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
#'   existing normalization can be \emph{updated} for a different unit of
#'   expression or after a conversion to a related spectral quantity. The
#'   accepted arguments for parameter \code{norm} are:
#'
#'   \describe{
#'   \item{\code{"max"}}{Scale all spectral values so that the value at the tallest peak is 1.0.}
#'   \item{\code{"min"}}{Scale all spectral values so that the value at the deepest valley is 1.0.}
#'   \item{\code{<number>}}{Scale all spectral values so that the value at \code{w.length == <number>} is 1.0, interpolating nearest wavelengths.}
#'   \item{\code{"update"}}{Normalise the data reusing the normalization criteria stored as metadata in attribute \code{normalization}.}
#'   \item{\code{"undo"}}{Reverse the normalisation of \code{x} restoring the original values of the spectral variable using the normalization metadata stored in attribute \code{normalization}.}
#'   \item{\code{"skip"}}{Return \code{x} unmodified.}
#'   }
#'
#'   In 'photobiology' (>= 0.10.8) detailed information about the normalization
#'   is stored in attribute \code{normalization}. \strong{In the case objects
#'   normalized using 'photobiology' (< 0.10.8), \code{norm} actions \code{"update"} and
#'   \code{"undo"} are not supported.}
#'
#'   In 'photobiology' (>= 0.10.10) applying a new normalization to an already
#'   normalized spectrum recomputes the multiplier factors stored in attribute
#'   \code{normalization} if it present. This ensures that the data in the
#'   returned object is the original one when \emph{undoing} multiple previous
#'   normalizations. Loss of precision due to floating-point arithmetic is
#'   possible, as well as in the format of attributes.
#'
#'   By default \code{max()} or \code{min()} are applied to the whole spectrum,
#'   but by passing a range of wavelengths as argument to parameter
#'   \code{tange}, they are applied to within this range of wavelengths.
#'   However, in all cases the normalization multipliers are applied to the
#'   whole spectrum.
#'
#'   By default the argument passed to \code{x} contains one or more \code{NA}
#'   values and the normalization is based on a summary quantity, the returned
#'   spectrum will contain only \code{NA} values. If \code{na.rm == TRUE} then
#'   the summary quantity will be calculated after striping \code{NA} values,
#'   and only the values that were \code{NA} in \code{x} will be \code{NA}
#'   values in the returned spectrum.
#'
#'   When the spectrum passed as argument to \code{x} had been previously
#'   scaled, in 'photobiology' (<= 0.10.9) the scaling attribute was always
#'   removed and no normalization factors returned. In 'photobiology'
#'   (>= 0.10.10) scaling information can be preserved by passing
#'   \code{keep.scaling = TRUE}.
#'
#'   When a numeric value is passed as argument to keep.scaling, the scaling
#'   uses \code{f = "total"} or \code{f = "mean"} depending on the class of
#'   \code{x}. Rescaling is only occasionally needed.
#'
#'   Method \code{normalize} is implemented for \code{solute_spct} objects but
#'   as the spectral data stored in them are a description of an intensive
#'   property of a substance, normalization is unlikely to useful. To represent
#'   solutions of specific concentrations of solutes, \code{filter_spct} objects
#'   should be used instead.
#'
#' @note The second formal argument is ellipsis, thus all parameters except
#'   \code{x} have to be always passed by name.
#'
#' @param x An R object
#' @param ... not used in current version
#'
#' @return A copy of the object passed as argument to \code{x} with the values
#'   of the spectral quantity rescaled to 1 at the normalization wavelength. If
#'   the normalization wavelength is not already present in \code{x}, it is
#'   added by interpolation---i.e. the returned value may be one row longer than
#'   \code{x} for each spectrum. Attributes \code{normalized} and
#'   \code{normalization} are set to keep a log of the computations applied.
#'   Other attributes are preserved.
#'
#' @examples
#' norm_sun.spct <- normalize(sun.spct, norm = "max")
#' str(is_normalized(norm_sun.spct))
#' str(getNormalization(norm_sun.spct))
#'
#' norm_sun.spct <- normalize(sun.spct, norm = 500)
#' str(is_normalized(norm_sun.spct))
#' str(getNormalization(norm_sun.spct))
#'
#' @export
#'
#' @family rescaling functions
#'
normalize <- function(x, ...) UseMethod("normalize")

#' @rdname normalize
#'
#' @note \code{normalise()} is another name for \code{normalize()}.
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
#'   length 2 with the limits of a range of wavelengths in nm. See
#'   \strong{Details}.
#' @param norm numeric Normalization wavelength (nm) or character string
#'   \code{"max"}, or \code{"min"} for normalization at the corresponding
#'   wavelength, \code{"update"} to update the normalization with previous
#'   criterion, \code{"undo"} to revert an existing normalization or
#'   \code{"skip"} to force return of \code{x} unchanged. See \strong{Details}.
#' @param unit.out No longer supported and is ignored with a warning.
#' @param keep.scaling logical or numeric Flag to indicate if any existing
#'   scaling should be preserved or not. The default, \code{FALSE}, preserves
#'   the behaviour of 'photobiology' (<= 0.10.9). See \strong{Details}.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before calculating the summary (e.g. \code{max()}) used for normalization.
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
#' @param col.names character vector containing the names of columns of \code{x}
#'   to be normalized. Other columns are retained unchanged.
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
#' @param .parallel	if \code{TRUE}, apply function in parallel, using parallel
#'   backend provided by packege 'foreach'.
#' @param .paropts a list of additional options passed into the \code{foreach}
#'   function when parallel computation is enabled. This is important if (for
#'   example) your code relies on external data or packages: use the
#'   \code{.export} and \code{.packages} arguments to supply them so that all
#'   cluster nodes have the correct environment set up for computing.
#'
#' @export
#'
#' @examples
#' norm_sun_evening.mspct <- normalize(sun_evening.mspct[1:3])
#' str(is_normalized(norm_sun_evening.mspct))
#' str(getNormalization(norm_sun_evening.mspct))
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
#' @note While method \code{getNormalized()} returns in some cases numeric
#'   values or possibly even character values if stored in attribute
#'   \code{normalized}, \code{is_normalized()} always returns a logical value
#'   and can be safely used in conditional code clauses.
#'
#' @export
#'
#' @family rescaling functions
#'
#' @examples
#' norm_sun.spct <- normalize(sun.spct)
#' is_normalized(norm_sun.spct)
#'
#' norm_sun_evening.mspct <- normalize(sun_evening.mspct[1:3])
#' str(is_normalized(norm_sun_evening.mspct))
#'
is_normalized <- function(x) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    spct.attr <- attr(x, "normalized", exact = TRUE)
    # in some versions a logical was used, but later the normalization wavelength
    # in old versions the attribute was set only when normalization was applied
    # for spectra in long form the attribute value is a named list.
    stopifnot(is.null(spct.attr) || is.numeric(spct.attr) ||
                is.logical(spct.attr) || is.list(spct.attr))
    return(!is.null(spct.attr) && any(as.logical(spct.attr)))
  } else if (is.generic_mspct(x)) {
    return(mslply(x, is_normalized))
  } else if (is.waveband(x)) {
    return(!is.na(normalization(x)))
  } else {
    warning("Method 'is_normalized()' not implemented for class: '",
            class(x)[1], "'.")
    return(NA)
  }
}

#' @rdname is_normalized
#'
#' @note \code{is_normalised()} is another name for \code{is_normalized()}.
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
#' @details In some computations spectral data that have been normalized need to
#'   be handled differently than data expressed in original units. Method
#'   \code{getNormalized()} makes it possible to query if a
#'   \code{generic_spct} or \code{generic_mspct} object or objects of derived
#'   classes contain data expressed in true physical units or normalized.
#'
#'   Method \code{getNormalization()} retrieves from objects storing spectral
#'   data, metadata that trace previously applied normalizations, making it
#'   possible to revert the effect of earlier normalizations. The metadata are
#'   also used when printing and plotting the spectra.
#'
#' @return For objects containing a single spectrum, \code{getNormalized()}
#'   returns a logical value, and exceptionally for objects created with
#'   'photobiology' (< 0.10.8), a numeric value (normalization wavelength
#'   expressed in nanometres). A character value for objects created with
#'   'photobiology' (< 0.x.x).  For collections of spectra and multiple spectra
#'   in long form, a named list, with one logical member for each spectrum is
#'   returned. If \code{x} is not a \code{generic_spct} object, \code{NA} is
#'   returned with a warning.
#'
#'   For objects containing a single spectrum, \code{getNormalization()} returns
#'   a list with five fields: \code{norm.type}, \code{norm.wl},
#'   \code{norm.factors}, \code{norm.cols}, \code{norm.range}. For collections
#'   of spectra, a named list of lists, with one member list for each member of
#'   the collection of spectra is returned. Objects created with versions of
#'   package 'photobiology' (< 0.10.8) are lacking the detailed normalization
#'   metadata, in which case \code{getNormalized()} must be used to detect
#'   normalization. See \code{\link{setNormalized}()} for the values stored in
#'   the fields. If \code{x} is not a \code{generic_spct} object, a named list
#'   with all fields set to \code{NA} is returned with a warning.
#'
#' @note While method \code{getNormalized()} returns in some cases numeric
#'   values or possibly even character values if stored in attribute
#'   \code{normalized}, \code{is_normalized()} always returns a logical value
#'   and can be safely used in conditional code clauses.
#'
#' @export
#'
#' @examples
#'
#' getNormalized(sun.spct)
#' str(getNormalization(sun.spct))
#'
#' norm_sun.spct <- normalize(sun.spct)
#'
#' is_normalized(norm_sun.spct)
#' getNormalized(norm_sun.spct)
#' str(getNormalization(norm_sun.spct))
#'
#' str(getNormalization(e2q(norm_sun.spct)))
#'
#' norm_gel.spct <- normalize(yellow_gel.spct)
#'
#' is_normalized(norm_gel.spct)
#' getNormalized(norm_gel.spct)
#' str(getNormalization(norm_gel.spct))
#'
#' getNormalization(T2Afr(norm_gel.spct))
#' getNormalization(any2A(norm_gel.spct))
#'
#' norm_sun_evening.mspct <- normalize(sun_evening.mspct[1:3])
#' str(is_normalized(norm_sun_evening.mspct))
#' str(getNormalized(norm_sun_evening.mspct))
#' str(getNormalization(norm_sun_evening.mspct))
#'
#' @family rescaling functions
#'
getNormalized <- function(x,
                          .force.numeric = FALSE) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    normalized <- attr(x, "normalized", exact = TRUE)
    if (is.null(normalized) || all(is.na(normalized))) {
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
  if (is.logical(normalized) && all(normalized)) {
    # backwards compatibility
    normalization.attr <- attr(x, "normalization")
    if (length(normalization.attr)) {
      norm.wl <- normalization.attr[["norm.wl"]]
      if (length(norm.wl) == 1L && is.numeric(norm.wl) && is.finite(norm.wl)) {
        normalized <- norm.wl
      }
    }
  }
  if (.force.numeric) {
    normalized <- suppressWarnings(as.numeric(normalized))
    if (!length(normalized)) {
      NA_real_
    } else {
      normalized
    }
  } else {
    normalized
  }
}

#' @rdname getNormalized
#'
#' @note \code{getNormalised()} is another name for \code{getNormalized()}.
#'
#' @export
#'
getNormalised <- getNormalized

#' @rdname getNormalized
#'
#' @export
#'
getNormalization <- function(x) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    if (is_normalized(x)) {
      # attribute in use >= 0.10.8
      normalization.list <- attr(x, "normalization", exact = TRUE)
      if (is.list(normalization.list)) {
        if (getMultipleWl(x) == 1L) {
          # check validity
          if (!exists("norm.range", normalization.list) &&
              exists("norm.type", normalization.list)) {
            # norm.range is missing 0.10.8 and 0.10.9
            normalization.list[["norm.range"]] <- rep(NA_real_, 2)
          }
          if (!exists("norm.wl", normalization.list) ||
              !is.numeric(normalization.list[["norm.wl"]])) {
            # norm.wl missing or corrupted breaks 'ggspectra'
            normalization.list[["norm.wl"]] <- NA_real_
          }
          # when removing columns the normalization data can remain behind
          # we return normalization only for existing columns
          selector <- normalization.list[["norm.cols"]] %in% colnames(x)
          if (!is.null(normalization.list[["norm.cols"]]) &&
              !any(is.na(normalization.list[["norm.cols"]])) &&
              !all(selector)) {
            normalization.list[["norm.type"]] <-
              normalization.list[["norm.type"]][selector]
            normalization.list[["norm.wl"]] <-
              normalization.list[["norm.wl"]][selector]
            normalization.list[["norm.factors"]] <-
              normalization.list[["norm.factors"]][selector]
            normalization.list[["norm.cols"]] <-
              normalization.list[["norm.cols"]][selector]
          }
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
#' @note \code{getNormalisation()} is another name for \code{getNormalization()}.
#'
getNormalisation <- getNormalization

#' Set the "normalized" and "normalization" attributes
#'
#' Function to write the \code{"normalized"} attribute of an existing
#' \code{generic_spct} object.
#'
#' @param x a generic_spct object.
#' @param norm logical or numeric If \code{FALSE} or \code{0} existing
#'   normalization metadata is is deleted from \code{x}. Otherwise, a positive
#'   numeric value is interpreted as the normalization wavelength (nm,
#'   nanometres).
#' @param norm.type character Type of normalization applied.
#' @param norm.factors numeric The scaling factor(s) so that dividing the
#'   spectral values by this factor reverts the normalization.
#' @param norm.cols character The name(s) of the columns in \code{x} that have
#'   been normalized.
#' @param norm.range numeric The wavelength range used for normalization (nm).
#' @param verbose logical Flag enabling or silencing informative warnings.
#'
#' @details This function \strong{is used internally}, although occasionally
#'   users may want to use it to "pretend" that spectral data have not been
#'   normalized. Use \code{\link{normalize}()} methods to apply a normalization
#'   and simultaneously set the metadata attributes. Function
#'   \code{setNormalized()} only saves to the attributes the metadata
#'   corresponding to an already applied normalization. The metadata provides a
#'   \emph{trace} of the transformations applied to spectral data that makes it
#'   possible to \emph{undo} the normalization. The metadata is also used in
#'   other by functions in package 'ggspectra' when automatically generating
#'   axis labels.
#'
#'   If \code{norm = FALSE} is passed in the call any normalization metadata
#'   present in \code{x} are reset and, thus, \code{x} marked as not normalized,
#'   without undoing the effect of the normalization.
#'
#'   If \code{x} is not a \code{generic_spct} object, \code{x} is not modified.
#'
#' @return The object passed as argument to \code{x} is modified by reference by
#'   adding attributes \code{normalized} and \code{normalization}, unless the
#'   argument passed to \code{x} is an anonymous expression. A copy of the
#'   modified \code{x} is always returned invisibly, even when setting by
#'   reference fails.
#'
#'   With a single spectrum in \code{x}, attribute \code{normalization} is set
#'   to a named list, and in the case of multiple spectra in long form, it is
#'   set a named list or named lists, unless \code{norm = FALSE} is passed, in
#'   which case the \code{normalization} attribute is deleted if already
#'   present. The named list for each spectrum contains the fields:
#'
#'   \describe{
#'   \item{\code{norm.type}}{\code{character} vector of length one, one of "max", "min", "wavelength".}
#'   \item{\code{norm.wl}}{\code{numeric}, normalization wavelength in nanometres.}
#'   \item{\code{norm.factors}}{\code{numeric}, multiplier constans used to scale the normalized spectral data.}
#'   \item{\code{norm.cols}}{\code{character}, the name of the columns of \code{x} that have been normalized.}
#'   \item{\code{norm.range}}{\code{numeric} vector of length 2, with min and maximum \code{w.length} values used to constrain the normalization.}
#'   }
#'
#'   Spectral objects, to allow choice in the trade-off between storage space
#'   and computation effort, can contain multiple columns of spectral data
#'   (e.g., \code{s.e.irrad} and \code{s.q.irard}), each one with its own
#'   normalization. Only one value is allowed for \code{norm.type}, but
#'   \code{norm.wl} and \code{norm.factors} are vectors of the same length as
#'   \code{norm.cols} as their numeric values depend on the unit/base of
#'   expression.
#'
#'   Attribute \code{normalized} is set to \code{TRUE} unless
#'   \code{norm = FALSE} is passed.
#'
#' @note Passing a \code{logical} as argument to \code{norm} is deprecated but
#'   accepted silently for backwards compatibility.
#'
#' @export
#'
#' @examples
#' norm_sun.spct <- normalize(sun.spct, norm = "max")
#' is_normalized(norm_sun.spct)
#' # rarely useful: pretend that the spectrum has not been normalized
#' pretended_sun.spct <- setNormalized(norm_sun.spct, norm = FALSE)
#' is_normalized(pretended_sun.spct)
#'
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

  if (is.generic_mspct(x)) {
    warning("To apply 'setNormalized()' to members of a collection ",
            "apply it with 'msmsply()'")
  }
  selector <- which(!is.na(norm))
  if (length(selector)) {
    if (!(is.numeric(norm) || is.logical(norm))) {
      stop("'norm' must be numeric or logical, but it is '", mode(norm), "'",
           call. = FALSE)
    }
    if (is.numeric(norm) && any(!is.finite(norm[selector]) | norm[selector] < 1)) {
      stop("If numeric, 'norm' must be finite and positive, not: '",
           paste(norm[selector], collapse = "', '"), "'",
           call. = FALSE)
    }
  }

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
    attr(x, "normalized") <- TRUE
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
#' @note \code{setNormalised()} is another name for \code{setNormalized()}.
#'
#' @export
#'
setNormalised <- setNormalized

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
  if (is.null(norm) || all(is.na(norm))) {
    norm <- "skip"
  }
  norm <- unique(na.omit(norm))
  if (length(norm) != 1L) {
    warning("'norm' is not unique! Using '", norm[1], "' instead of '",
            paste(norm, collapse = "', '"), "'")
    norm <- norm[1]
  }
  if (!(is.numeric(norm) || is.character(norm))) {
    stop("'norm' must be numeric or character, not: '", mode(norm), "'",
         call. = FALSE)
  }
  if (is.numeric(norm) && any(!is.finite(norm) | norm < 1)) {
    stop("If numeric, 'norm' must be finite and positive, not: '", norm, "'",
         call. = FALSE)
  }
  valid.norm.values <-
    c("skip", "undo", "update", "max", "maximum", "min", "minimum", "wipe.attrs")
  if (is.character(norm) && !norm  %in% valid.norm.values) {
    stop("If character, 'norm' must be one of '",
         paste(valid.norm.values, collapse = "', '"),
         "', not: '", norm, "'",
         call. = FALSE)
  }

  scale.is.dirty <- FALSE

  # if 'norm' is a character vector, we use the first element
  # thus, all columns always get the same type of normalization
  if (is.character(norm) && length(norm) > 1) {
    if (length(unique(norm)) > 1) {
      warning("Multiple 'norm' values supplied by name. Using the first one: ",
              norm[1], ".")
    }
    norm <- norm[1]
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

  updating <- all(unlist(is_normalized(spct), use.names = FALSE))

  if (updating) {
    # we retrieve the existing normalization data
    old.normalization.ls <- getNormalization(spct)

    required.fields <- c("norm.type", "norm.wl", "norm.cols", "norm.range")
    has.normalization.metadata <-
      length(old.normalization.ls) >= length(required.fields) &&
      all(required.fields %in% names(old.normalization.ls)) &&
      !any(is.na(unlist(old.normalization.ls[required.fields])))

    if (!has.normalization.metadata && norm[1] %in% c("update", "undo")) {
      warning("Normalization not updated/undone: action not supported for ",
              "objects lacking normalization metadata.")
      return(spct)
    } else if (has.normalization.metadata) {
      if (norm[1] == "update") {
        # extract old normalization criteria
        if (old.normalization.ls$norm.type[1] == "wavelength") {
          norm <- old.normalization.ls$norm.wl
        } else {
          norm <- old.normalization.ls$norm.type
        }
        range <- old.normalization.ls$norm.range
      }

      # remove the old normalization
      if (norm[1] == "wipe.attrs") {
        # remove attributes
        spct <- denormalize_spct(spct,
                                 wipe.away = TRUE)
        return(spct)
      } else {
        # restore to original scale and remove attributes
        spct <- denormalize_spct(spct,
                                 wipe.away = FALSE)
        if (norm[1] == "undo") {
          return(spct)
        }
      }
    } else {
      spct <- denormalize_spct(spct,
                               wipe.away = TRUE)
      scale.is.dirty <- TRUE # no normalization factors stored in attribute
    }
  } else if (norm[1] %in% c("update", "undo", "wipe.attrs")) {
    # not normalized, nothing to update
    return(spct)
  }

  # Here spct contains spectral data reverted to original scaling

  # we do not remove NA's from the returned object, only from working copy x
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

  norm.arg <- norm # for later use

  # normalization will wipe out any existing scaling except for its effect
  # on the computed factors.
  if (is_scaled(x) && !keep.scaling) {
    # The behaviour in <= 0.10.9
    # remove scaling metadata and do not save norm.factors
    scale.is.dirty <- TRUE
    setScaled(spct, scaled = FALSE)
  } else {
    # retain scaling metadata and save norm.factors
    scale.is.dirty <- scale.is.dirty || FALSE
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
        # target normalization wavelength is within range
        norm.wl <- norm[i]
        tmp.spct <- spct[ , c("w.length", col)]
        class(tmp.spct) <- class(spct)
        scale.factor <- 1 /
          interpolate_spct(spct = tmp.spct, w.length.out = norm.wl)[ , eval(col)]
      } else {
        warning("'norm = ", norm[i], "' value(s) outside range of ",
                round(range[1], 1), " to ", round(range[2], 1), " (nm)")
        scale.factor <- NA_real_
        norm.wl <- NA_real_
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
  if (scale.is.dirty) {
    message("'norm.factors' not stored")
  }
  z # setNormalized makes its returned value invisible
}

#' Undo previously applied normalization
#'
#' @param spct A \code{generic_spct} object or a \code{generic_mspct} object,
#'   or objects of derived classes.
#' @param wipe.away logical If \code{TRUE} the normalization metadata is removed
#'   without undoing the effect of the normalization, and otherwise after
#'   restoring the spectral data to its original scaling.
#'
#' @return A modified copy of \code{spct} if it was previously normalized or
#' \code{spct} unchanged, otherwise.
#'
#' @keywords internal
#'
denormalize_spct <- function(spct, wipe.away = FALSE) {
  if (!all(unlist(is_normalized(spct), use.names = FALSE))) {
    return(spct)
  }
  # collection of spectra
  if (is.generic_mspct(spct)) {
    return(
      msmsply(mspct = spct, .fun = denormalize_spct, wipe.away = wipe.away)
    )
  }
  # if wiping away single spectrum or long form spectra
  if (wipe.away) {
    message("Removing normalization metadata keeping normalization!")
    attr(spct, "normalized") <- FALSE
    attr(spct, "normalization") <- NULL
    return(spct)
  }
  # if undoing normalization, done spectrum by spectrum
  if (is.generic_spct(spct) && getMultipleWl(spct) > 1L) {
    mspct <- subset2mspct(spct)
    mspct <- msmsply(mspct = mspct,
                    .fun = denormalize_spct,
                    wipe.away = wipe.away)
    return(rbindspct(mspct, idfactor = getIdFactor(spct)))
  }

  # undo normalization of a single spectrum
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

#' Restore normalization
#'
#' After altering the columns present in a spectrum restore normalization
#' with the previously used criteria.
#'
#' @param x generic_spct One or more spectra stored in long form, with no
#'   normalization applied.
#' @param old.normalization.ls list A list describing the normalization criteria
#'   in the format returned by \code{getNormalization()}.
#'
#' @details The normalization criteria are extracted from
#' \code{old.normalization.ls} and applied. In the case of spectra in
#' long form, the normalization can be stored as a named list with
#' values for each spectrum as members.
#'
#' If \code{x} contains multiple spectra in long form, the same normalization
#' criterion is applied to all of them based on \code{old.normalization.ls}.
#' If \code{getNormalization()} is a named list and stored values of
#' \code{norm.type} are equal, it is used. If \code{norm.type == "wavelength"},
#' \code{norm.wl} is used if consistent. \code{norm.range} is used as is if
#' consistent, and otherwise expanded to the joint range.
#' If metadata are no consistent across multiple spectra stored in long form
#' re-normalization is skipped with a warning.
#'
#' @keywords internal
#'
restore_normalization <- function(x, old.normalization.ls) {

  if (is_normalized(x)) {
    warning("'x' is already normalized. Skipping!!")
    return(x)
  }

  # apply the pre-existing normalization criteria
  # to columns present in changed x
  if (getMultipleWl(x) > 1L) {
    old.norm <- unique(as.vector(sapply(old.normalization.ls,
                                        `[[`,
                                        i = "norm.type",
                                        USE.NAMES = FALSE)))
    if (length(old.norm) > 1L) {
      warning("Inconsistent normalization type across ",
              getMultipleWl(x), " spectra. Skipping!!")
      return(x)
    }
  } else {
    old.norm <- old.normalization.ls$norm.type
  }

  if (old.norm[1] == "wavelength") {
    if (getMultipleWl(x) > 1L) {
      old.norm <- unique(as.vector(sapply(old.normalization.ls,
                                          `[[`,
                                          i = "norm.wl",
                                          USE.NAMES = FALSE)))
      if (length(old.norm) > 1L) {
        warning("Inconsistent normalization wavelength across ",
                getMultipleWl(x), " spectra. Skipping!!")
        return(x)
      }
    } else {
      old.norm <- old.normalization.ls$norm.wl
    }
  }
  if (getMultipleWl(x) > 1L) {
    old.mins <-
      as.vector(sapply(old.normalization.ls,
                       function(x) { x[["norm.range"]][1L] },
                       USE.NAMES = FALSE))
    old.maxs <-
      as.vector(sapply(old.normalization.ls,
                       function(x) { x[["norm.range"]][2L] },
                       USE.NAMES = FALSE))
    # use tolerance of 1 nm
    if (diff(range(old.mins)) > 1 || diff(range(old.maxs)) > 1) {
      warning("Inconsistent normalization range across ",
              getMultipleWl(x), " spectra. Using joint range!!")
      old.range <- wl_range(x)
    } else {
      old.range <- c(min(old.mins), max(old.maxs))
    }
  } else {
    old.range <- old.normalization.ls$norm.range
  }

  normalize(x,
            range = old.range,
            norm = old.norm,
            keep.scaling = TRUE)
}
