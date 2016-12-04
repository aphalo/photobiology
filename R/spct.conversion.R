#' Conversion from counts per second to physical quantities
#'
#' Conversion of spectral data expressed as cps into irradiance, transmittance
#' or reflectance.
#'
#' @param x.sample,x.clear,x.opaque,x.white,x.black cps_spct objects.
#' @param pre.fun function A function applied to x.sample before converison.
#' @param dyn.range numeric The effective dynamic range of the instrument,
#'    if \code{NULL} it is automatically set based on integartion time
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
  setSourceSpct(z)
  if (length(descriptor[["inst.calib"]][["wl.range"]]) == 2) {
    clip_wl(z, descriptor[["inst.calib"]][["wl.range"]])
  }
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
  stopifnot(getInstrDesc(x.sample)$spectrometer.sn ==
              getInstrDesc(x.white)$spectrometer.sn)

  instr.desc <- list(getInstrDesc(x.sample),
                     getInstrDesc(x.white))
  if (!is.null(x.black)) {
    instr.desc <- c(instr.desc, getInstrDesc(x.black))
  }

  if (anyNA(instr.desc)) {
    warning("Missing intrument descriptor attributes.")
  } else {
    instr.sn <- sapply(instr.desc, `[[`, i = "spectrometer.sn")
    if (!length(unique(instr.sn)) == 1) {
      stop("ERROR: serial number mismatch between cps_spct objects")
    }
  }

  cps.col.sample <- grep("^cps", names(x.sample), value = TRUE)
  cps.col.white <- grep("^cps", names(x.white), value = TRUE)
  stopifnot(length(cps.col.sample) == 1 && length(cps.col.white) == 1)
  other.cols <- setdiff(names(x.sample), cps.col.sample)
  z <- as.generic_spct(x.sample)
  z[[cps.col.sample]] <- NULL
  z[["Rfr"]] <- x.sample[[cps.col.sample]] / x.white[[cps.col.white]]
  # guess of dynamic range as a function of bracketing for sample
  if (is.null(dyn.range)) {
    acq_settings <- getInstrSettings(x.sample)
    if (is.na(acq_settings)) {
      dyn.range <- 7e2
    } else {
      integ.range <- range(acq_settings[["integ.time"]])
      dyn.range <- min(7e2 * integ.range[2] / integ.range[1], 1e4)
    }
  }
  # based on dynamic range and spectrum of light source we set bad data to NA
  z[["Rfr"]] <- ifelse(x.white[[cps.col.white]] <
                         (max(x.white[[cps.col.white]], na.rm = TRUE)) /
                         dyn.range,
                       NA_real_,
                       z[["Rfr"]])
  setReflectorSpct(z)
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
              is.cps_spct(x.clear) &&
               (is.null(x.opaque) || is.cps_spct(x.opaque)))
  instr.desc <- list(getInstrDesc(x.sample),
                 getInstrDesc(x.clear))
  if (!is.null(x.opaque)) {
    instr.desc <- c(instr.desc, getInstrDesc(x.opaque))
  }

  if (anyNA(instr.desc)) {
    warning("Missing intrument descriptor attributes.")
  } else {
    instr.sn <- sapply(instr.desc, `[[`, i = "spectrometer.sn")
    if (!length(unique(instr.sn)) == 1) {
      stop("ERROR: serial number mismatch between cps_spct objects")
    }
  }

  if (!is.null(x.opaque)) {
    x.sample <- x.sample - x.opaque
    x.clear <- x.clear - x.opaque
  }

  cps.col.sample <- grep("^cps", names(x.sample), value = TRUE)
  cps.col.clear <- grep("^cps", names(x.clear), value = TRUE)
  stopifnot(length(cps.col.sample) == 1 && length(cps.col.clear) == 1)
  z <- as.generic_spct(x.sample)
  z[[cps.col.sample]] <- NULL
  z[["Tfr"]] <- x.sample[[cps.col.sample]] / x.clear[[cps.col.clear]]

  # guess of dynamic range as a function of bracketing for sample
  if (is.null(dyn.range)) {
    acq_settings <- getInstrSettings(x.sample)
    if (is.na(acq_settings)) {
      dyn.range <- 7e2
    } else {
      integ.range <- range(acq_settings[["integ.time"]])
      dyn.range <- min(7e2 * integ.range[2] / integ.range[1], 1e4)
    }
  }
  # based on dynamic range and spectrum of light source we set bad data to NA
  z[["Tfr"]] <- ifelse(x.clear[[cps.col.clear]] <
                         (max(x.clear[[cps.col.clear]], na.rm = TRUE)) /
                         dyn.range,
                       NA_real_,
                       z[["Tfr"]])
  setFilterSpct(z)
}
