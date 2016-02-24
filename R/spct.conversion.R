#' Conversion from counts per second to physical quantities
#'
#' Conversion of spectral data expressed as cps into irradiance, transmittance
#' or reflectance.
#'
#' @param x.sample,x.clear,x.opaque,x.white,x.black cps_spct objects.
#' @param pre.fun function A function applied to x.sample before converison.
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
  irrad.mult <- getInstrDesc(x.sample)$inst.calib$irrad.mult
  if (!is.null(pre.fun)) {
    x.sample <- pre.fun(x.sample, ...)
  }
  cps.col <- grep("^cps", names(x.sample), value = TRUE)
  stopifnot(length(cps.col) == 1)
  other.cols <- setdiff(names(x.sample), cps.col)
  z <- as.generic_spct(x.sample)[other.cols]
  z[["s.e.irrad"]] <- x.sample[[cps.col]] * irrad.mult
  z[[cps.col]] <- NULL
  setSourceSpct(z)
}

#' @rdname cps2irrad
#' @export
cps2Rfr <- function(x.sample, x.white, x.black = NULL) {
  # we make sure that all input spectra have been measured with the same
  # instrument by comparing serial numbers
  stopifnot(is.cps_spct(x.sample) &&
              !is.null(getInstrDesc(x.sample)))
  stopifnot(is.cps_spct(x.white) &&
              !is.null(getInstrDesc(x.white)))
  stopifnot(getInstrDesc(x.sample)$spectrometer.sn ==
              getInstrDesc(x.white)$spectrometer.sn)
  if (!is.null(x.black)) {
    stopifnot(is.cps_spct(x.black) &&
                !is.null(getInstrDesc(x.black)))
    stopifnot(getInstrDesc(x.sample)$spectrometer.sn ==
                getInstrDesc(x.black)$spectrometer.sn)
    x.sample <- x.sample - x.black
    x.white <- x.white - x.black
  }
  cps.col.sample <- grep("^cps", names(x.sample), value = TRUE)
  cps.col.white <- grep("^cps", names(x.white), value = TRUE)
  stopifnot(length(cps.col.sample) == 1 && length(cps.col.white) == 1)
  other.cols <- setdiff(names(x.sample), cps.col.sample)
  z <- as.generic_spct(x.sample)[other.cols]
  z[["Rfr"]] <- x.sample[[cps.col.sample]] / x.white[[cps.col.white]]
  z[[cps.col.sample]] <- NULL
  setReflectorSpct(z)
}

#' @rdname cps2irrad
#' @export
cps2Tfr <- function(x.sample, x.clear, x.opaque = NULL) {
  # we make sure that all input spectra have been measured with the same
  # instrument by comparing serial numbers
  stopifnot(is.cps_spct(x.sample) &&
              !is.null(getInstrDesc(x.sample)))
  stopifnot(is.cps_spct(x.clear) &&
              !is.null(getInstrDesc(x.clear)))
  stopifnot(getInstrDesc(x.sample)$spectrometer.sn ==
              getInstrDesc(x.clear)$spectrometer.sn)
  if (!is.null(x.opaque)) {
    stopifnot(is.cps_spct(x.opaque) &&
                !is.null(getInstrDesc(x.opaque)))
    stopifnot(getInstrDesc(x.sample)$spectrometer.sn ==
                getInstrDesc(x.opaque)$spectrometer.sn)
    x.sample <- x.sample - x.opaque
    x.clear <- x.clear - x.opaque
  }
  cps.col.sample <- grep("^cps", names(x.sample), value = TRUE)
  cps.col.clear <- grep("^cps", names(x.clear), value = TRUE)
  stopifnot(length(cps.col.sample) == 1 && length(cps.col.clear) == 1)
  other.cols.sample <- setdiff(names(x.sample), cps.col.sample)
  z <- as.generic_spct(x.sample)[other.cols.sample]
  z[["Rfr"]] <- x.sample[[cps.col.sample]] / x.clear[[cps.col.clear]]
  z[["Rfr"]] <- ifelse(x.clear[[cps.col.clear]] < 1e-3 * max(x.clear[[cps.col.clear]]),
                       NA_real_,
                       z[["Rfr"]])
  z[[cps.col.sample]] <- NULL
  setReflectorSpct(z)
}
