#' Conversion from counts per second to physical quantities
#'
#' Conversion of spectral data expressed as cps into irradiance, transmittance
#' or reflectance.
#'
#' @param x.sample,x.clear,x.opaque,x.white,x.black cps_spct objects.
#' @param pre.fun function A function applied to x.sample before conversion.
#' @param dyn.range numeric The effective dynamic range of the instrument,
#'    if \code{NULL} it is automatically set based on integration time
#'    bracketing.
#' @param ... Additional arguments passed to \code{pre.fun}.
#'
#' @return A source_spct, filter_spct or reflector_spct object containing the
#'   spectral values expressed in physical units.
#'
#' @note In contrast to other classes defined in package 'photobiology', class
#'   "cps_spct" can have more than one column of cps counts in cases where the
#'   intention is to merge these values as part of the processing at the time
#'   the calibration is applied. However, being these functions the final step
#'   in the conversion to physical units, they accept as input only objects
#'   with a single "cps" column, as merging is expected to have been already
#'   done.
#'
#' @export
#'
cps2irrad <- function(x.sample, pre.fun = NULL, ...) {
  stopifnot(is.cps_spct(x.sample) &&
              !is.null(getInstrDesc(x.sample)) &&
              !is.null(getInstrSettings(x.sample)))
  descriptor <- getInstrDesc(x.sample)
  irrad.mult <- descriptor[["inst.calib"]][["irrad.mult"]]
  if (!is.null(pre.fun)) {
    x.sample <- pre.fun(x.sample, ...)
  }
  cps.col.sample <- grep("^cps", names(x.sample), value = TRUE)
  stopifnot(length(cps.col.sample) == 1)
  z <- as.generic_spct(x.sample)
  z[[cps.col.sample]] <- NULL
  z[["s.e.irrad"]] <- x.sample[[cps.col.sample]] * irrad.mult
  setSourceSpct(z, time.unit = getTimeUnit(x.sample))
  z <- copy_attributes(x.sample, z)
  if (length(descriptor[["inst.calib"]][["wl.range"]]) == 2) {
    z <- clip_wl(z, descriptor[["inst.calib"]][["wl.range"]])
  }
  z
}

#' @rdname cps2irrad
#' @export
cps2Rfr <- function(x.sample,
                    x.white,
                    x.black = NULL,
                    dyn.range = NULL) {
  # we make sure that all input spectra have been measured with the same
  # instrument by comparing serial numbers
  stopifnot(is.cps_spct(x.sample) &&
              !is.null(getInstrDesc(x.sample)))
  stopifnot(is.cps_spct(x.white) &&
              !is.null(getInstrDesc(x.white)))
  stopifnot(getInstrDesc(x.sample)[["spectrometer.sn"]] ==
            getInstrDesc(x.white)[["spectrometer.sn"]])

  if (!is.null(x.black)) {
    stopifnot(getInstrDesc(x.sample)[["spectrometer.sn"]] ==
                getInstrDesc(x.black)[["spectrometer.sn"]])
  }

#  instr.desc <- getInstrDesc(x.sample)

  cps.col.sample <- grep("^cps", names(x.sample), value = TRUE)
  cps.col.white <- grep("^cps", names(x.white), value = TRUE)
  stopifnot(length(cps.col.sample) == 1L && length(cps.col.white) == 1)
#  other.cols <- setdiff(names(x.sample), cps.col.sample)
  z <- as.generic_spct(x.sample)
  z[[cps.col.sample]] <- NULL
  if (!is.null(x.black)) {
    cps.col.black <- grep("^cps", names(x.black), value = TRUE)
    stopifnot(length(cps.col.black) == 1L)
    z[["Rfr"]] <- (x.sample[[cps.col.sample]] - x.black[[cps.col.black]])  /
      (x.white[[cps.col.white]] - x.black[[cps.col.black]])
  } else {
    z[["Rfr"]] <- x.sample[[cps.col.sample]] / x.white[[cps.col.white]]
  }
  # guess of dynamic range as a function of bracketing for sample
  if (is.null(dyn.range)) {
    acq_settings <- getInstrSettings(x.sample)
    if (!is.list(acq_settings) && is.na(acq_settings)) {
      dyn.range <- 7e2
    } else {
      integ.time <- acq_settings[["integ.time"]]
      dyn.range <- min(7e2 * max(integ.time) / min(integ.time), 1e4)
      dyn.range <- dyn.range * acq_settings[["rel.signal"]]
    }
  }
  # based on dynamic range and spectrum of light source we set bad data to NA
  z[["Rfr"]] <- ifelse(x.white[[cps.col.white]] <
                         (max(x.white[[cps.col.white]], na.rm = TRUE)) /
                         dyn.range,
                       NA_real_,
                       z[["Rfr"]])
  setReflectorSpct(z)
  copy_attributes(x.sample, z)
}

#' @rdname cps2irrad
#' @export
cps2Tfr <- function(x.sample,
                    x.clear,
                    x.opaque = NULL,
                    dyn.range = NULL) {
  # we make sure that all input spectra have been measured with the same
  # instrument by comparing serial numbers
  stopifnot(is.cps_spct(x.sample) &&
              !is.null(getInstrDesc(x.sample)))
  stopifnot(is.cps_spct(x.clear) &&
              !is.null(getInstrDesc(x.clear)))
  stopifnot(getInstrDesc(x.sample)[["spectrometer.sn"]] ==
              getInstrDesc(x.clear)[["spectrometer.sn"]])

  if (!is.null(x.opaque)) {
    stopifnot(getInstrDesc(x.sample)[["spectrometer.sn"]] ==
                getInstrDesc(x.opaque)[["spectrometer.sn"]])
  }

#  instr.desc <- getInstrDesc(x.sample)

  cps.col.sample <- grep("^cps", names(x.sample), value = TRUE)
  cps.col.clear <- grep("^cps", names(x.clear), value = TRUE)
  stopifnot(length(cps.col.sample) == 1 && length(cps.col.clear) == 1)
  z <- as.generic_spct(x.sample)
  z[[cps.col.sample]] <- NULL
  if (!is.null(x.opaque)) {
    cps.col.opaque <- grep("^cps", names(x.opaque), value = TRUE)
    stopifnot(length(cps.col.opaque) == 1L)
    z[["Tfr"]] <- (x.sample[[cps.col.sample]] - x.opaque[[cps.col.opaque]])  /
      (x.clear[[cps.col.clear]] - x.opaque[[cps.col.opaque]])
  } else {
    z[["Tfr"]] <- x.sample[[cps.col.sample]] / x.clear[[cps.col.clear]]
  }

  # guess of dynamic range as a function of bracketing for sample
  if (is.null(dyn.range)) {
    acq_settings <- getInstrSettings(x.sample)
    if (!is.list(acq_settings) && is.na(acq_settings)) {
      dyn.range <- 7e2
    } else {
      integ.time <- acq_settings[["integ.time"]]
      dyn.range <- min(7e2 * max(integ.time) / min(integ.time), 1e4)
      dyn.range <- dyn.range * acq_settings[["rel.signal"]]
    }
  } else {
    warning("Spectrometer settings are not available.\n",
            "Pass a suitable value as argument to 'dyn.range'.\n")
  }
  # based on dynamic range and spectrum of light source we set bad data to NA
  z[["Tfr"]] <- ifelse(x.clear[[cps.col.clear]] <
                         (max(x.clear[[cps.col.clear]], na.rm = TRUE)) /
                         dyn.range,
                       NA_real_,
                       z[["Tfr"]])
  setFilterSpct(z)
  copy_attributes(x.sample, z)
}
