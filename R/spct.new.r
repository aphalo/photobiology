# Constructors ------------------------------------------------------------

#' Spectral-object constructors
#'
#' These constructor functions can be used to create spectral objects derived
#' from \code{generic_spct}. They take as arguments numeric vectors for the
#' wavelengths and spectral data, and numeric, character, and logical values for
#' metadata attributes to be saved to the objects created and options
#' controlling the creation process.
#'
#' @param w.length numeric vector with wavelengths in nanometres [\eqn{nm}].
#' @param s.e.irrad numeric vector with spectral energy irradiance in
#'   [\eqn{W\,m^{-2}\,nm^{-1}}] or [\eqn{J\,d^{-1}\,m^{-2}\,nm^{-1}}{J d-1 m-2 nm-1}].
#' @param s.q.irrad numeric A vector with spectral photon irradiance in
#'   [\eqn{mol\,s^{-1}\,m^{-2}\,nm^{-1}}{mol s-1 m-2 nm-1}] or
#'   [\eqn{mol\,d^{-1}\,m^{-2}\,nm^{-1}}{mol d-1 m-2 nm-1}].
#' @param time.unit character string indicating the time unit used for spectral
#'   irradiance or exposure (\code{"second"}, \code{"day"} or \code{"exposure"})
#'   or an object of class duration as defined in package lubridate.
#' @param bswf.used character A string indicating the BSWF used, if any, for
#'   spectral effective irradiance or exposure (\code{"none"} or the name of the
#'   BSWF).
#' @param comment character A string to be added as a comment attribute to the
#'   object created.
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning.
#' @param multiple.wl	numeric Maximum number of repeated \code{w.length} entries
#'   with same value. (As with multiple spectra stored in long from).
#' @param idfactor character Name of factor distinguishing multiple spectra when
#'   stored longitudinally (required if \code{multiple.wl} > 1).
#' @param ... other arguments passed to \code{tibble()} such as vectors or
#'   factors to be added as additional columns.
#'
#' @return A object of class \code{generic_spct} or a class derived from it,
#'   depending on the function used. In other words an object of a class with
#'   the same name as the constructor function.
#'
#' @details Constructors can be used to create spectral objects from spectral
#'   quantities expressed on a single base or unit. Some of the functions have
#'   different formal parameters accepting a quantity expressed in different
#'   units, however, an argument can be passed to only one of these formal
#'   parameters in a given call. The constructors \code{object_spct()} and
#'   \code{chroma_spct()} require arguments to be passed for multiple but
#'   distinct spectral quantities.
#'
#' @export
#'
#' @family constructors of spectral objects
#'
#' @rdname source_spct
#'
source_spct <- function(w.length = NULL,
                        s.e.irrad = NULL,
                        s.q.irrad = NULL,
                        ...,
                        time.unit = c("second", "day", "exposure"),
                        bswf.used = c("none", "unknown"),
                        comment = NULL,
                        strict.range = getOption("photobiology.strict.range", default = FALSE),
                        multiple.wl = 1L,
                        idfactor = NULL) {
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
#' @param irrad.mult numeric vector with multipliers for each detector pixel
#'   expressed in units of \eqn{W\,m^{-2}\,nm^{-1}\,n^{-1}\,s}{W m-2 nm-1 n-1 s},
#'   where \eqn{n\,s^{-1}}{n s-1} are detector counts per second.
#'
#' @export
#'
calibration_spct <- function(w.length = NULL,
                             irrad.mult = NA_real_,
                             ...,
                             comment = NULL,
                             instr.desc = NA,
                             multiple.wl = 1L,
                             idfactor = NULL) {
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
#' @param counts numeric vector with raw counts expressed per scan.
#' @param instr.desc a list describing the spectrometer used to acquire the data.
#' @param instr.settings a list describing the settings used to acquire the data.
#'
#' @export
#'
raw_spct <- function(w.length = NULL,
                     counts = NA_real_,
                     ...,
                     comment = NULL,
                     instr.desc = NA,
                     instr.settings = NA,
                     multiple.wl = 1L,
                     idfactor = NULL) {
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
#'   [\eqn{n\,s^{-1}}{n s-1}]
#'
#' @export
#'
cps_spct <- function(w.length = NULL,
                     cps = NA_real_,
                     ...,
                     comment = NULL,
                     instr.desc = NA,
                     instr.settings = NA,
                     multiple.wl = 1L,
                     idfactor = NULL) {
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
                         ...,
                         comment = NULL,
                         multiple.wl = 1L,
                         idfactor = NULL) {
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
#' @param s.e.response numeric vector with a biological, chemical or physical
#'   response expressed per unit spectral energy irradiance
#'   [\eqn{W\,m^{-2}\,nm^{-1}}{W m-2 nm-1} or \eqn{J\,d^{-1}\,m^{-2}\,nm^{-1}}{J d-1 m-2 nm-1}].
#' @param s.q.response numeric vector with a biological, chemical or physical
#'   response expressed per unit spectral photon irradiance in
#'   [\eqn{mol\,s^{-1}\,m^{-2}\,nm^{-1}}{mol s-1 m-2 nm-1} or \eqn{mol\,d^{-1}\,m^{-2}\,nm^{-1}}{mol d-1 m-2 nm-1}].
#' @param response.type a character string, either \code{"response"} or
#'   \code{"action"}.
#'
#' @export
#'
response_spct <- function(w.length = NULL,
                          s.e.response = NULL,
                          s.q.response = NULL,
                          ...,
                          time.unit = c("second", "day", "exposure"),
                          response.type = c("response", "action"),
                          comment = NULL,
                          multiple.wl = 1L,
                          idfactor = NULL) {
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
                  time.unit = time.unit,
                  response.type = response.type,
                  multiple.wl = multiple.wl,
                  idfactor = idfactor)
  z
}

#' @rdname source_spct
#'
#' @param Tfr numeric vector with spectral transmittance as fraction of one
#'   [\eqn{/1}].
#' @param Tpc numeric vector with spectral transmittance as percent values
#' @param Afr numeric vector of absorptance as fraction of one [\eqn{/1}].
#' @param A   numeric vector of absorbance values (\eqn{log_{10}}{log10}-base
#'   a.u.)
#' @param Tfr.type character string indicating whether transmittance and
#'   absorptance values are \code{"total"} or \code{"internal"} values
#' @param Rfr.constant numeric The value of the reflection factor [\eqn{/1}].
#' @param thickness numeric The thickness of the material.
#' @param attenuation.mode character One of \code{"reflection"},
#'   \code{"absorption"} or \code{"mixed"}.
#'
#' @section Warning for filter_spct!: Not entering metadata when creating an
#'   object will limit the available operations! While "internal" transmittance
#'   is defined as the transmittance of the material body itself, "total"
#'   transmittance includes the effects of surface reflectance on the amount of
#'   light transmitted. For non-diffusing materials like glass an approximate
#'   \code{Rfr.constant} value can be used to convert "total" into "internal"
#'   transmittance values and vice versa. Use \code{NA} if not known, or not
#'   applicable, e.g., for materials subject to internal scattering.
#'
#' @seealso \code{\link{setFilterProperties}}
#'
#' @export
#'
filter_spct <- function(w.length = NULL,
                        Tfr = NULL,
                        Tpc = NULL,
                        Afr = NULL,
                        A = NULL,
                        ...,
                        Tfr.type = c("total", "internal"),
                        Rfr.constant = NA_real_,
                        thickness = NA_real_,
                        attenuation.mode = NA,
                        comment = NULL,
                        strict.range = getOption("photobiology.strict.range", default = FALSE),
                        multiple.wl = 1L,
                        idfactor = NULL) {
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
                Rfr.constant = Rfr.constant,
                thickness = thickness,
                attenuation.mode = attenuation.mode,
                strict.range = strict.range,
                multiple.wl = multiple.wl,
                idfactor = idfactor)
  z
}

#' @rdname source_spct
#'
#' @param Rfr numeric vector with spectral reflectance as fraction of one
#'   [\eqn{/1}].
#' @param Rpc numeric vector with spectral reflectance as percent values.
#' @param Rfr.type character A string, either \code{"total"} or
#'   \code{"specular"}.
#'
#' @export
#'
reflector_spct <- function(w.length = NULL,
                           Rfr = NULL,
                           Rpc = NULL,
                           ...,
                           Rfr.type = c("total", "specular"),
                           comment = NULL,
                           strict.range = getOption("photobiology.strict.range", default = FALSE),
                           multiple.wl = 1L,
                           idfactor = NULL) {
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
#' @param K.mole numeric vector with molar attenuation coefficient in SI units
#'    [\eqn{m^2\,mol^-1}{m2 mol-1}].
#' @param K.mass numeric vector with mass attenuation coefficient in SI units
#'    [\eqn{m^2\,g^-1}{m2 g-1}].
#' @param attenuation.XS numeric vector with attenuation cross section values
#'    (Converted during object construction into \code{K.mole}.)
#' @param K.type character A string, either \code{"attenuation"},
#'   \code{"absorption"} or \code{"scattering"}.
#' @param mass numeric The molar mass in Dalton [Da] (\eqn{Da = g\,mol^{-1}}{Da = g mol-1}).
#' @param formula character The molecular formula.
#' @param structure raster A bitmap of the structure.
#' @param name,solvent.name character The names of the substance and of the
#'   solvent. A named character vector, with member names such as "IUPAC" for
#'   the authority.
#' @param ID,solvent.ID character The ID of the substance and of the solvent. A
#'   named character vector, with member names such as "ChemSpider" or "PubChem"
#'   for the authority.
#' @param log.base numeric Normally one of \code{e} or \code{10}. Data are
#'   stored always on base 10 corresponding to decadal absorbance as used in
#'   chemistry.
#'
#' @section Warning for solute_spct!:
#'   You should always set the base for logarithms to match that on which the
#'   absorbance data are expressed. Failing to do this will result in bad data
#'   and all further computation will be wrong. Not entering metadata when
#'   creating an object will limit the available operations! Mass should be
#'   indicated in daltons or \eqn{g\,mol^{-1}}{g mol-1}. The SI unit of molar attenuation
#'   coefficient is the square metre per mole (\eqn{m^2\,mol^{1}}{m2 mol-1}),
#'   but in practice, quantities are usually expressed in terms of
#'   \eqn{M^{-1}\,cm^{-1}} or \eqn{l\,mol^{-1}\,cm^{-1}} (the latter two units are
#'   both equal to 0.1 \eqn{m^2\,mol^{-1}} and quantities expressed in them need
#'   to be divided by 10 when passed as arguments to \code{K.mole}.).
#'
#' @seealso \code{\link{setSoluteProperties}}
#'
#' @export
#'
solute_spct <- function(w.length = NULL,
                        K.mole = NULL,
                        K.mass = NULL,
                        attenuation.XS = NULL,
                        ...,
                        log.base = 10,
                        K.type = c("attenuation", "absorption", "scattering"),
                        name = NA_character_,
                        mass = NA_character_,
                        formula = NULL,
                        structure = grDevices::as.raster(matrix()),
                        ID = NA_character_,
                        solvent.name = NA_character_,
                        solvent.ID = NA_character_,
                        comment = NULL,
                        strict.range = getOption("photobiology.strict.range", default = FALSE),
                        multiple.wl = 1L,
                        idfactor = NULL) {
  if (!is.null(attenuation.XS) && is.null(K.mole)) {
    K.mole <- attenuation.XS / 3.82343216e-21 # epsilon = sigma * N_A / (log(10) * 1e3)
  }
  if (length(w.length) == 0) {
    z <- tibble::tibble(w.length = numeric(), K.mole = numeric(), ...)
  } else if (is.null(K.mass) && is.numeric(K.mole)) {
    if (log.base != 10) {
      K.mole <- log10(K.mole^log.base)
    }
    z <- tibble::tibble(w.length, K.mole, ...)
  } else if (is.null(K.mole) && is.numeric(K.mass)) {
    stop("Support for 'K.mass' not yet implemented.")
    if (log.base != 10) {
      K.mass <- log10(K.mass^log.base)
    }
    z <- tibble::tibble(w.length, K.mass, ...)
  } else {
    warning("Only one of K.mole, or K.mass should be different from NULL.")
    z <- tibble::tibble(w.length, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setSoluteSpct(x = z,
                K.type = K.type,
                name = name,
                mass = mass,
                formula = formula,
                structure = structure,
                ID = ID,
                solvent.name = solvent.name,
                solvent.ID = solvent.ID,
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
                        ...,
                        Tfr.type = c("total", "internal"),
                        Rfr.type = c("total", "specular"),
                        comment = NULL,
                        strict.range = getOption("photobiology.strict.range", default = FALSE),
                        multiple.wl = 1L,
                        idfactor = NULL) {
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
                        ...,
                        comment = NULL,
                        strict.range = getOption("photobiology.strict.range", default = FALSE),
                        multiple.wl = 1L,
                        idfactor = NULL) {
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
#' @param x an R object.
#' @param ... other arguments passed to "set" functions.
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
#' @param x an R object.
#' @param ... other arguments passed to "set" functions.
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
#' @param x an R object.
#' @param ... other arguments passed to "set" functions.
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
#' @param x an R object.
#' @param time.unit character string indicating the time unit used for spectral
#'   irradiance or exposure (\code{"second"}, \code{"day"} or \code{"exposure"})
#'   or an object of class duration as defined in package lubridate.
#' @param bswf.used character A string indicating the BSWF used, if any, for
#'   spectral effective irradiance or exposure (\code{"none"} or the name of the
#'   BSWF).
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning.
#' @param ... other arguments passed to "set" functions.
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
           time.unit = c("second", "day", "exposure"),
           bswf.used = c("none", "unknown"),
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
#' @param x an R object.
#' @param time.unit character string indicating the time unit used for spectral
#'   irradiance or exposure (\code{"second"}, \code{"day"} or \code{"exposure"})
#'   or an object of class duration as defined in package lubridate.
#' @param ... other arguments passed to "set" functions.
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
#' @param x an R object.
#' @param Rfr.type a character string, either \code{"total"} or
#'   \code{"specular"}.
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning.
#' @param ... other arguments passed to "set" functions.
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
#' @param x an R object.
#' @param Tfr.type a character string, either \code{"total"} or
#'   \code{"internal"}.
#' @param Rfr.type a character string, either \code{"total"} or
#'   \code{"specular"}.
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning.
#' @param ... other arguments passed to "set" functions.
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

#' Coerce or convert into a filter spectrum
#'
#' Return a possibly modified copy of an R object with its class set to a filter
#' spectrum. In the case of conversion from a \code{solute_spct} object, compute
#' the spectral quantity based on additional input from user.
#'
#' @param x an R object.
#' @param Tfr.type a character string, either \code{"total"} or
#'   \code{"internal"}.
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning.
#' @param ... other arguments passed to "set" functions.
#'
#' @return A copy of \code{x} converted into a \code{filter_spct}. object.
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

#'@describeIn as.filter_spct
#'
#' @param Rfr.constant numeric The value of the reflection factor (/1) to be
#'   set.
#' @param comment character A string to be added as a comment attribute to the
#'   object created. If not supplied, the comment will be copied from \code{x}.
#' @param molar.concentration,mass.concentration numeric Concentration to be
#'   used to compute transmittance of the solute in solution
#'   [\eqn{mol\,m^{-3} = mmol\,dm^{-3}}{mol m-3 = mmol dm-3} or
#'   \eqn{kg\,m^{-3} = g\,dm^{-3}}{kg m-3 = g dm-3}, respectively].
#' @param path.length numeric The length of the light path (\eqn{m}) used to
#'   compute transmittance of the solute in a solution.
#'
#' @export
#'
as.filter_spct.solute_spct <-
  function(x,
           Tfr.type = "internal",
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           Rfr.constant = NA_real_,
           comment = NULL,
           molar.concentration = NULL,
           mass.concentration = NULL,
           path.length = 1, # meter
           ...) {
    stopifnot(Tfr.type == "internal" ||
                Tfr.type == "total" && !is.na(Rfr.constant) && Rfr.constant >= 0)
    stopifnot(xor(is.null(molar.concentration), is.null(mass.concentration)))
    solute.properties <- getSoluteProperties(x)
    attenuation.mode <- ifelse(getKType(x) == "absorption",
                               "absorption",
                               "mixed")
    # we do calculations using moles
    if (is.null(molar.concentration)) {
      molar.concentration <- mass.concentration / solute.properties[["mass"]]
    }
    if (! "K.mole" %in% colnames(x)) {
      x[["K.mole"]] <- x[["K.mass"]] / solute.properties[["mass"]]
    }
    if (is.null(comment)) {
      comment <- paste("Computed from 'solute_spct' for ",
                       solute.properties[["name"]], ".\n",
                       comment(x), sep = "")
    }
    z <- filter_spct(w.length = x[["w.length"]],
                    A = x[["K.mole"]] * molar.concentration * path.length,
                    Tfr.type = "internal",
                    Rfr.constant = Rfr.constant,
                    thickness = path.length,
                    attenuation.mode = attenuation.mode,
                    comment = comment,
                    strict.range = strict.range,
                    multiple.wl = getMultipleWl(x),
                    ...)
    other.cols <-
      setdiff(colnames(x), c("w.length", "K.mole", "K.mass"))
    if (length(other.cols)) {
      zz <- cbind(z, x[ , other.cols])
      copy_attributes(z, zz)
    } else {
      z
    }
  }

#' Coerce to a solute spectrum
#'
#' Return a possibly modified copy of an R object with its class set to
#' \code{solute_spct} (a solute spectrum). In the case of conversion from a
#' \code{filter_spct} object, compute spectral molar attenuation  based on
#' additional input from user.
#'
#' @param x an R object.
#' @param K.type a character string, one of \code{"attenuation"},
#'   \code{"absorption"} or \code{"scattering"}.
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning.
#' @param ... other arguments passed to "set" functions.
#'
#' @return A copy of \code{x} converted into a \code{solute_spct} object.
#'
#' @seealso \code{\link{setSoluteSpct}}
#'
#' @export
#'
#' @family constructors of spectral objects
#'
as.solute_spct <- function(x, ...) {UseMethod("as.solute_spct")}

#' @describeIn as.solute_spct
#'
#'
#' @export
#'
as.solute_spct.default <-
  function(x,
           K.type = c("attenuation", "absorption", "scattering"),
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           ...) {
    setSoluteSpct(x,
                  K.type = K.type,
                  strict.range = strict.range,
                  ...)
  }

#'@describeIn as.solute_spct
#'
#' @param name,solvent.name character The names of the substance and of the
#'   solvent. A named character vector, with member names such as "IUPAC" for
#'   the authority.
#' @param mass numeric The mass in Dalton (Da = g/mol).
#' @param formula character The molecular formula.
#' @param structure raster A bitmap of the structure.
#' @param ID,solvent.ID character The IDs of the substance and of the solvent. A
#'   named character vector, with member names such as "ChemSpider" or "PubChen"
#'   for the authority.
#' @param comment character A string to be added as a comment attribute to the
#'   object created. If not supplied, the comment will be copied from \code{x}.
#' @param molar.concentration,mass.concentration numeric Concentration to be
#'   used to compute transmittance of the solute in solution [\eqn{mol\,m^{-3} =
#'   mmol\,dm^{-3}}{mol m-3 = mmol dm-3} or \eqn{kg\,m^{-3} = g\,dm^{-3}}{kg m-3
#'   = g dm-3}, respectively].
#' @param path.length numeric The length of the light path (\eqn{m}) used to
#'   compute transmittance of the solute in a solution.
#'
#' @export
#'
as.solute_spct.filter_spct <-
  function(x,
           K.type = c("attenuation", "absorption", "scattering"),
           name = NA_character_,
           mass = NA_character_,
           formula = NULL,
           structure = grDevices::as.raster(matrix()),
           ID = NA_character_,
           solvent.name = NA_character_,
           solvent.ID = NA_character_,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           comment = NULL,
           molar.concentration = NULL,
           mass.concentration = NULL,
           path.length = 1, # meter
           ...) {
    if (getTfrType(x) != "internal") {
      x <- convertTfrType(x, Tfr.type = "internal")
    }
    x <- any2A(x, action = "replace")
    # we do calculations using moles
    if (is.null(molar.concentration)) {
      molar.concentration <- mass.concentration / mass
    }
    if (is.null(comment)) {
      comment <- paste("Computed from 'filter_spct' for ",
                       name, ".\n",
                       comment(x), sep = "")
    }
    z <- solute_spct(w.length = x[["w.length"]],
                     K.mole = x[["A"]] / molar.concentration / path.length,
                     log.base = 10,
                     K.type = K.type,
                     name = name,
                     mass = mass,
                     formula = formula,
                     structure = structure,
                     comment = comment,
                     ID = ID,
                     solvent.name = solvent.name,
                     solvent.ID = solvent.ID,
                     strict.range = strict.range,
                     multiple.wl = getMultipleWl(x),
                     ...)
    other.cols <-
      setdiff(colnames(x), c("w.length", "A"))
    if (length(other.cols)) {
      zz <- cbind(z, x[ , other.cols])
      copy_attributes(z, zz)
    } else {
      z
    }
  }

#' Coerce to a spectrum
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object.
#' @param ... other arguments passed to "set" functions.
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
#' Merge a \code{filter_spct} with a \code{reflector_spct} returning an
#' \code{object_spct} object, even if wavelength values are mismatched.
#'
#' @param x,y a \code{filter_spct} object and a \code{reflector_spct} object.
#' @param by a vector of shared column names in \code{x} and \code{y} to merge
#'   on; \code{by} defaults to \code{w.length}.
#' @param ... other arguments passed to \code{dplyr::inner_join()}.
#' @param w.length.out numeric vector of wavelengths to be used for the returned
#'   object (\eqn{nm}).
#' @param Tfr.type.out character string indicating whether transmittance values
#'   in the returned object should be expressed as \code{"total"} or
#'   \code{"internal"}. This applies only to the case when an \code{object_spct}
#'   is returned.
#'
#' @note If a numeric vector is supplied as argument for \code{w.length.out},
#'   the two spectra are interpolated to the new wavelength values before
#'   merging. The default argument for \code{w.length.out} is
#'   \code{x[["w.length"]]}.
#'
#' @return An \code{object_spct} is returned as the result of merging a
#'   \code{filter_spct} and a \code{reflector_spct} object.
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
