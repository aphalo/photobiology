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
  stopifnot(is.cps_spct(x.sample))
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
  if (!is.null(x.black)) {
    x.sample <- x.sample - x.black
    x.white <- x.white - x.black
  }
  stopifnot(is.cps_spct(x.sample) &&
              is.cps_spct(x.white))
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
  if (!is.null(x.opaque)) {
    x.sample <- x.sample - x.opaque
    x.clear <- x.clear - x.opaque
  }
  stopifnot(is.cps_spct(x.sample) &&
              is.cps_spct(x.clear))
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