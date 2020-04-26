

# Constructors ------------------------------------------------------------

#' Spectral-object constructor
#'
#' These functions can be used to create spectral objects derived from
#' \code{generic_spct}. They take as arguments numeric vectors for the data
#' character scalars for attributes, and a logical flag.
#'
#' @param w.length numeric vector with wavelengths in nanometres
#' @param s.e.irrad numeric vector with spectral energy irradiance in [W m-2
#'   nm-1] or [J d-1 m-2 nm-1]
#' @param s.q.irrad numeric A vector with spectral photon irradiance in [mol s-1
#'   m-2 nm-1] or [mol d-1 m-2 nm-1].
#' @param time.unit character string indicating the time unit used for spectral
#'   irradiance or exposure ("second" , "day" or "exposure") or an object of
#'   class duration as defined in package lubridate.
#' @param bswf.used character A string indicating the BSWF used, if any, for
#'   spectral effective irradiance or exposure ("none" or the name of the BSWF).
#' @param comment character A string to be added as a comment attribute to the
#'   object created.
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning.
#' @param multiple.wl	numeric Maximum number of repeated w.length entries with
#'   same value.
#' @param idfactor character Name of factor distinguishing multiple spectra when
#'   stored logitudinally (required if mulitple.wl > 1).
#' @param ... other arguments passed to \code{tibble()}
#'
#' @return A object of class generic_spct or a class derived from it, depending
#'   on the function used. In other words an object of a class with the same
#'   name as the constructor function.
#'
#' @export
#'
#' @family constructors of spectral objects
#'
#' @rdname source_spct
#'
#' @note The functions can be used to add only one spectral quantity to a
#'   spectral object. Some of the functions have different arguments, for the
#'   same quantity expressed in different units. An actual parameter can be
#'   supplied to only one of these formal parameters in a given call to any of
#'   these functions.
#'
source_spct <- function(w.length = NULL,
                        s.e.irrad = NULL,
                        s.q.irrad = NULL,
                        time.unit = c("second", "day", "exposure"),
                        bswf.used = c("none", "unknown"),
                        comment = NULL,
                        strict.range = getOption("photobiology.strict.range", default = FALSE),
                        multiple.wl = 1L,
                        idfactor = NULL,
                        ...) {
  if (length(w.length) == 0) {
    z <- tibble::tibble(w.length = numeric(), s.e.irrad = numeric(), ...)
  } else if (is.null(s.q.irrad) && (is.numeric(s.e.irrad))) {
    z <- tibble::tibble(w.length, s.e.irrad, ...)
  } else if (is.null(s.e.irrad) && (is.numeric(s.q.irrad))) {
    z <- tibble::tibble(w.length, s.q.irrad, ...)
  } else {
    warning("Only one of s.e.irrad or s.q.irrad should be different from NULL.")
    z <- tibble::tibble(w.length, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setSourceSpct(z,
                time.unit = time.unit,
                bswf.used = bswf.used,
                strict.range = strict.range,
                multiple.wl = multiple.wl,
                idfactor = idfactor)
  z
}

#' @rdname source_spct
#'
#' @param irrad.mult numeric vector with multipliers for each detector pixel.
#'
#' @export
#'
calibration_spct <- function(w.length = NULL,
                             irrad.mult = NA_real_,
                             comment = NULL,
                             instr.desc = NA,
                             multiple.wl = 1L,
                             idfactor = NULL,
                             ...) {
  if (length(w.length) == 0) {
    z <- tibble::tibble(w.length = numeric(), irrad.mult = numeric(), ...)
  } else {
    z <- tibble::tibble(w.length = w.length, irrad.mult = irrad.mult, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setCalibrationSpct(z,
                     multiple.wl = multiple.wl,
                     idfactor = idfactor)
  setInstrDesc(z, instr.desc)
  z
}

#' @rdname source_spct
#'
#' @param counts numeric vector with raw counts expressed per scan
#' @param instr.desc a list
#' @param instr.settings a list
#'
#' @export
#'
raw_spct <- function(w.length = NULL,
                     counts = NA_real_,
                     comment = NULL,
                     instr.desc = NA,
                     instr.settings = NA,
                     multiple.wl = 1L,
                     idfactor = NULL,
                     ...) {
  if (length(w.length) == 0) {
    z <- tibble::tibble(w.length = numeric(), counts = numeric(), ...)
  } else {
    z <- tibble::tibble(w.length = w.length, counts = counts, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setRawSpct(z,
             multiple.wl = multiple.wl,
             idfactor = idfactor)
  setInstrDesc(z, instr.desc)
  setInstrSettings(z, instr.settings)
  z
}

#' @rdname source_spct
#'
#' @param cps numeric vector with linearized raw counts expressed per second
#'
#' @export
#'
cps_spct <- function(w.length = NULL,
                     cps = NA_real_,
                     comment = NULL,
                     instr.desc = NA,
                     instr.settings = NA,
                     multiple.wl = 1L,
                     idfactor = NULL,
                     ...) {
  if (any(grepl("^cps", names(list(...)))) && is.na(cps)) {
    if (length(w.length) == 0) {
      z <- tibble::tibble(w.length = numeric(), ...)
    } else {
      z <- tibble::tibble(w.length = w.length, ...)
    }
  } else {
    if (length(w.length) == 0) {
      z <- tibble::tibble(w.length = numeric(), cps = numeric(), ...)
    } else {
      z <- tibble::tibble(w.length = w.length, cps = cps, ...)
    }
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setCpsSpct(z,
             multiple.wl = multiple.wl,
             idfactor = idfactor)
  setInstrDesc(z, instr.desc)
  setInstrSettings(z, instr.settings)
  z
}

#' @rdname source_spct
#'
#' @export
#'
generic_spct <- function(w.length = NULL,
                         comment = NULL,
                         multiple.wl = 1L,
                         idfactor = NULL,
                         ...) {
  if (length(w.length) == 0) {
    z <- tibble::tibble(w.length = numeric(), ...)
  } else {
    z <- tibble::tibble(w.length = w.length, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setGenericSpct(z,
                 multiple.wl = multiple.wl,
                 idfactor = idfactor)
  z
}

#' @rdname source_spct
#'
#' @param s.e.response numeric vector with spectral energy irradiance in W m-2
#'   nm-1 or J d-1 m-2 nm-1
#' @param s.q.response numeric vector with spectral photon irradiance in mol s-1
#'   m-2 nm-1 or mol d-1 m-2 nm-1
#'
#' @export
#'
response_spct <- function(w.length = NULL,
                          s.e.response = NULL,
                          s.q.response = NULL,
                          time.unit = c("second", "day", "exposure"),
                          comment = NULL,
                          multiple.wl = 1L,
                          idfactor = NULL,
                          ...) {
  if (length(w.length) == 0) {
    z <- tibble::tibble(w.length = numeric(), s.e.response = numeric(), ...)
  } else if (is.null(s.q.response) && (is.numeric(s.e.response))) {
    z <- tibble::tibble(w.length, s.e.response, ...)
  } else if (is.null(s.e.response) && (is.numeric(s.q.response))) {
    z <- tibble::tibble(w.length, s.q.response, ...)
  } else {
    warning("Only one of s.e.response or s.q.response should be different from NULL.")
    z <- tibble::tibble(w.length, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setResponseSpct(z,
                  time.unit,
                  multiple.wl = multiple.wl,
                  idfactor = idfactor)
  z
}

#' @rdname source_spct
#'
#' @param Tfr numeric vector with spectral transmittance as fraction of one
#' @param Tpc numeric vector with spectral transmittance as percent values
#' @param Afr numeric vector of absorptance as fraction of one
#' @param A   numeric vector of absorbance values (log10 based a.u.)
#' @param Tfr.type character string indicating whether transmittance
#'   and absorptance values are "total" or "internal" values
#'
#' @note "internal" transmittance is defined as the transmittance of the
#'   material body itself, while "total" transmittance includes the effects of
#'   surface reflectance on the amount of light transmitted.
#'
#' @export
#'
filter_spct <- function(w.length = NULL,
                        Tfr = NULL,
                        Tpc = NULL,
                        Afr = NULL,
                        A = NULL,
                        Tfr.type = c("total", "internal"),
                        comment = NULL,
                        strict.range = getOption("photobiology.strict.range", default = FALSE),
                        multiple.wl = 1L,
                        idfactor = NULL,
                        ...) {
  if (length(w.length) == 0) {
    z <- tibble::tibble(w.length = numeric(), Tfr = numeric())
  } else if (is.null(Tpc) && is.null(A) && is.null(Afr) && is.numeric(Tfr)) {
    z <- tibble::tibble(w.length, Tfr, ...)
  } else if (is.null(Tfr) && is.null(A) && is.null(Afr) && is.numeric(Tpc)) {
    z <- tibble::tibble(w.length, Tpc, ...)
  } else if (is.null(Tpc) && is.null(Tfr) && is.null(Afr) && is.numeric(A)) {
    z <- tibble::tibble(w.length, A, ...)
  } else if (is.null(Tpc) && is.null(Tfr) && is.null(A) && is.numeric(Afr)) {
    z <- tibble::tibble(w.length, Afr, ...)
  } else {
    warning("Only one of Tfr, Tpc, Afr, or A should be different from NULL.")
    z <- tibble::tibble(w.length, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setFilterSpct(x = z,
                Tfr.type = Tfr.type,
                strict.range = strict.range,
                multiple.wl = multiple.wl,
                idfactor = idfactor)
  z
}

#' @rdname source_spct
#'
#' @param Rfr numeric vector with spectral reflectance as fraction of one
#' @param Rpc numeric vector with spectral reflectance as percent values
#' @param Rfr.type character A string, either "total" or "specular".
#'
#' @export
#'
reflector_spct <- function(w.length = NULL,
                           Rfr=NULL,
                           Rpc=NULL,
                           Rfr.type=c("total", "specular"),
                           comment=NULL,
                           strict.range = getOption("photobiology.strict.range", default = FALSE),
                           multiple.wl = 1L,
                           idfactor = NULL,
                           ...) {
  if (length(w.length) == 0) {
    z <- tibble::tibble(w.length = numeric(), Rfr = numeric(), ...)
  } else if (is.null(Rpc) && is.numeric(Rfr)) {
    z <- tibble::tibble(w.length, Rfr, ...)
  } else if (is.null(Rfr) && is.numeric(Rpc)) {
    z <- tibble::tibble(w.length, Rpc, ...)
  } else {
    warning("Only one of Rfr, or Rpc should be different from NULL.")
    z <- tibble::tibble(w.length, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setReflectorSpct(x = z,
                   Rfr.type = Rfr.type,
                   strict.range = strict.range,
                   multiple.wl = multiple.wl,
                   idfactor = idfactor)
  z
}

#' @rdname source_spct
#'
#' @export
#'
object_spct <- function(w.length = NULL,
                        Rfr = NULL,
                        Tfr = NULL,
                        Afr = NULL,
                        Tfr.type = c("total", "internal"),
                        Rfr.type = c("total", "specular"),
                        comment = NULL,
                        strict.range = getOption("photobiology.strict.range", default = FALSE),
                        multiple.wl = 1L,
                        idfactor = NULL,
                        ...) {
  if (length(w.length) == 0) {
    z <- tibble::tibble(w.length = numeric(),
                        Rfr = numeric(), Tfr = numeric(), ...)
  } else if (is.null(Afr)) {
    z <- tibble::tibble(w.length, Rfr, Tfr, ...)
  } else if (is.null(Tfr)) {
    z <- tibble::tibble(w.length, Rfr, Afr, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setObjectSpct(z,
                Tfr.type = Tfr.type,
                Rfr.type = Rfr.type,
                strict.range = strict.range,
                multiple.wl = multiple.wl,
                idfactor = idfactor)
  z
}

#' @rdname source_spct
#'
#' @param x,y,z numeric colour coordinates
#'
#' @export
#'
chroma_spct <- function(w.length=NULL,
                        x,
                        y,
                        z,
                        comment=NULL,
                        strict.range = getOption("photobiology.strict.range", default = FALSE),
                        multiple.wl = 1L,
                        idfactor = NULL,
                        ...) {
  if (length(w.length) == 0) {
    z <- tibble::tibble(w.length = numeric(),
                        x = numeric(), y = numeric(), z = numeric(), ...)
  } else {
    z <- tibble::tibble(w.length, x, y, z, ...)
    if (!is.null(comment)) {
      comment(z) <- comment
    }
    setChromaSpct(z,
                  multiple.wl = multiple.wl,
                  idfactor = idfactor)
  }
  z
}

# as methods for spct classes --------------------------------------------

# defined as generics + default as additional methods are defined in other
# packages of the R for Photobiology suite.

#' Coerce to a spectrum
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object
#' @param ... other arguments passed to "set" functions
#'
#' @return A copy of \code{x} converted into a \code{generic_spct} object.
#'
#' @seealso \code{\link{setGenericSpct}}
#'
#' @export
#'
#' @family constructors of spectral objects
#'
as.generic_spct <- function(x, ...) {UseMethod("as.generic_spct")}

#' @describeIn as.generic_spct
#'
#' @export
#'
as.generic_spct.default <- function(x, ...) {
  setGenericSpct(x, ...)
}

#' Coerce to a spectrum
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object
#' @param ... other arguments passed to "set" functions
#'
#' @return A copy of \code{x} converted into a \code{calibration_spct} object.
#'
#' @seealso \code{\link{setGenericSpct}}
#'
#' @export
#'
#' @family constructors of spectral objects
#'
as.calibration_spct <- function(x, ...) {UseMethod("as.calibration_spct")}

#' @describeIn as.calibration_spct
#'
#' @export
#'
as.calibration_spct.default <- function(x, ...) {
  setCalibrationSpct(x, ...)
}

#' Coerce to a spectrum
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object
#' @param ... other arguments passed to "set" functions
#'
#' @return A copy of \code{x} converted into a \code{raw_spct} object.
#'
#' @seealso \code{\link{setGenericSpct}}
#'
#' @export
#'
#' @family constructors of spectral objects
#'
as.raw_spct <- function(x, ...) {UseMethod("as.raw_spct")}

#' @describeIn as.raw_spct
#'
#' @export
#'
as.raw_spct.default <- function(x, ...) {
  setRawSpct(x, ...)
}

#' Coerce to a spectrum
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object
#' @param ... other arguments passed to "set" functions
#'
#' @return A copy of \code{x} converted into a \code{cps_spct} object.
#'
#' @seealso \code{\link{setGenericSpct}}
#'
#' @export
#'
#' @family constructors of spectral objects
#'
as.cps_spct <- function(x, ...) {UseMethod("as.cps_spct")}

#' @describeIn as.cps_spct
#'
#' @export
#'
as.cps_spct.default <- function(x, ...) {
  setCpsSpct(x, ...)
}

#' Coerce to a spectrum
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object
#' @param time.unit character A string, "second", "day" or "exposure"
#' @param bswf.used character
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning
#' @param ... other arguments passed to "set" functions
#'
#' @return A copy of \code{x} converted into a \code{source_spct} object.
#'
#' @seealso \code{\link{setGenericSpct}}
#'
#' @export
#'
#' @family constructors of spectral objects
#'
as.source_spct <- function(x, ...) {UseMethod("as.source_spct")}

#' @describeIn as.source_spct
#'
#' @export
#'
as.source_spct.default <-
  function(x,
           time.unit=c("second", "day", "exposure"),
           bswf.used=c("none", "unknown"),
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           ...) {
    setSourceSpct(x,
                  time.unit = time.unit,
                  strict.range = strict.range,
                  bswf.used = bswf.used,
                  ...)
  }

#' Coerce to a spectrum
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object
#' @param time.unit character A string, "second", "day" or "exposure"
#' @param ... other arguments passed to "set" functions
#'
#' @return A copy of \code{x} converted into a \code{response_spct} object.
#'
#' @seealso \code{\link{setGenericSpct}}
#'
#' @export
#'
#' @family constructors of spectral objects
#'
as.response_spct <- function(x, ...) {UseMethod("as.response_spct")}

#' @describeIn as.response_spct
#'
#' @export
#'
as.response_spct.default <- function(x, time.unit = "second", ...) {
  setResponseSpct(x, time.unit = time.unit, ...)
}

#' Coerce to a spectrum
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object
#' @param Tfr.type a character string, either "total" or "internal"
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning
#' @param ... other arguments passed to "set" functions
#'
#' @return A copy of \code{x} converted into a \code{filter_spct} object.
#'
#' @seealso \code{\link{setGenericSpct}}
#'
#' @export
#'
#' @family constructors of spectral objects
#'
as.filter_spct <- function(x, ...) {UseMethod("as.filter_spct")}

#'@describeIn as.filter_spct
#'
#' @export
#'
as.filter_spct.default <-
  function(x,
           Tfr.type = c("total", "internal"),
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           ...) {
    setFilterSpct(x,
                  Tfr.type = Tfr.type,
                  strict.range = strict.range,
                  ...)
  }

#' Coerce to a spectrum
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object
#' @param Rfr.type a character string, either "total" or "specular"
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning
#' @param ... other arguments passed to "set" functions
#'
#' @return A copy of \code{x} converted into a \code{reflector_spct} object.
#'
#' @seealso \code{\link{setGenericSpct}}
#'
#' @export
#'
#' @family constructors of spectral objects
#'
as.reflector_spct <- function(x, ...) {UseMethod("as.reflector_spct")}

#' @describeIn as.reflector_spct
#'
#'
#' @export
#'
as.reflector_spct.default <-
  function(x,
           Rfr.type = c("total", "specular"),
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           ...) {
    setReflectorSpct(x,
                     Rfr.type = Rfr.type,
                     strict.range = strict.range,
                     ...)
  }

#' Coerce to a spectrum
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object
#' @param Tfr.type a character string, either "total" or "internal"
#' @param Rfr.type a character string, either "total" or "specular"
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning
#' @param ... other arguments passed to "set" functions
#'
#' @return A copy of \code{x} converted into a \code{object_spct} object.
#'
#' @seealso \code{\link{setGenericSpct}}
#'
#' @export
#'
#' @family constructors of spectral objects
#'
as.object_spct <- function(x, ...) {UseMethod("as.object_spct")}

#' @describeIn as.object_spct
#'
#' @export
#'
as.object_spct.default <- function(x,
                                   Tfr.type=c("total", "internal"),
                                   Rfr.type=c("total", "specular"),
                                   strict.range = getOption("photobiology.strict.range", default = FALSE),
                                   ...) {
  setObjectSpct(x,
                Tfr.type = Tfr.type,
                Rfr.type = Rfr.type,
                strict.range = strict.range,
                ...)
}

#' Coerce to a spectrum
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object
#' @param ... other arguments passed to "set" functions
#'
#' @return A copy of \code{x} converted into a \code{chroma_spct} object.
#'
#' @seealso \code{\link{setGenericSpct}}
#'
#' @export
#'
#' @family constructors of spectral objects
#'
as.chroma_spct <- function(x, ...) {UseMethod("as.chroma_spct")}

#' @describeIn as.chroma_spct
#'
#' @export
#'
as.chroma_spct.default <- function(x, ...) {
  setChromaSpct(x, ...)
}

# merge -------------------------------------------------------------------


#' Merge into object_spct
#'
#' Merge a filter_spct with a reflector_spct returning an object_spct object,
#' even if wavelength values are mismatched.
#'
#' @param x,y a filter_spct object and a reflector_spct object.
#' @param by a vector of shared column names in \code{x} and \code{y} to merge
#'   on; \code{by} defaults to \code{w.length}.
#' @param ... other arguments passed to \code{dplyr::inner_join()}
#' @param w.length.out numeric vector of wavelengths to be used for the returned
#'   object (nm).
#' @param Tfr.type.out character string indicating whether transmittance
#'   values in the returned object should be expressed as "total" or "internal".
#'   This applies only to the case when an object_spct is returned.
#'
#' @note If a numeric vector is supplied as argument for \code{w.length.out},
#'   the two spectra are interpolated to the new wavelength values before
#'   merging. The default argument for \code{w.length.out} is x[[w.length]].
#'
#' @return An object_spct is returned as the result of merging a filter_spct and
#'   a reflector_spct object.
#'
#' @seealso \code{\link[dplyr]{join}}
#'
#' @export
#'
merge2object_spct <- function(x, y,
                              by = "w.length", ...,
                              w.length.out =  x[["w.length"]],
                              Tfr.type.out = "total") {
  class.x <- class(x)
  class.y <- class(y)
  stopifnot(("filter_spct" %in% class.x && "reflector_spct" %in% class.y) ||
              ("reflector_spct" %in% class.x && "filter_spct" %in% class.y))
  stopifnot(!is.unsorted(w.length.out, strictly = TRUE))

  if ("filter_spct" %in% class.x && "reflector_spct" %in% class.y) {
    xx <- x
    yy <- y
  } else {
    xx <- y
    yy <- x
  }

  xx <- any2T(xx, action = "replace")
  xx <- interpolate_spct(spct = xx,
                         w.length.out = w.length.out,
                         fill = NA,
                         length.out = NULL)
  yy <- interpolate_spct(spct = yy,
                         w.length.out = w.length.out,
                         fill = NA,
                         length.out = NULL)

  # strip class attributes so that within dplyr code they are tibbles
  Tfr.type <- getTfrType(xx)
  Rfr.type <- getRfrType(yy)
  rmDerivedSpct(xx)
  rmDerivedSpct(yy)
  z <- dplyr::inner_join(xx, yy, by = "w.length", ...)
  if (Tfr.type.out == "internal" && Tfr.type == "total") {
    stopifnot(Rfr.type == "total")
    z[["Tfr"]] <- z[["Tfr"]] / (1 - z[["Rfr"]])
    Tfr.type <- "internal"
  }
  if (Tfr.type.out == "total" && Tfr.type == "internal") {
    stopifnot(Rfr.type == "total")
    z[["Tfr"]] <- z[["Tfr"]] * (1 - z[["Rfr"]])
    Tfr.type <- "total"
  }
  setObjectSpct(z, Tfr.type = Tfr.type, Rfr.type = Rfr.type)
  zz <- merge_attributes(x, y, z, which.not = c("Tfr.type", "Rfr.type", "comment"))

  comment(zz) <- paste("Merged spectrum\ncomment(x):\n",
                      comment(x),
                      "\nclass: ",
                      class_spct(x)[1L],
                      "\n\ncomment(y):\n",
                      comment(y),
                      "\nclass: ",
                      class_spct(y)[1L])
  zz
}
