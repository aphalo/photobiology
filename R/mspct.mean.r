#' mean
#'
#' Mean method for collections of spectra. Computes the mean at each wavelength
#' across all the spectra in the collection.
#'
#' @param x An R object. Currently this package defines methods for collections of
#'    spectral objects.
#' @param trim	numeric. The fraction (0 to 0.5) of observations to be trimmed from
#'   each end of x before the mean is computed. Values of trim outside that
#'   range are taken as the nearest endpoint.
#' @param na.rm	logical. A value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param ...	Further arguments passed to or from other methods.
#'
#' @note Trimming of extreme values and omission of NAs is done separately at
#' each wavelength. Interpolation is not applied, so all spectra in \code{x}
#' must share the same set of wavelengths.
#'
#' @seealso \code{\link[base]{mean}}
#'
#' @name mean
#'
#' @export
#'

mean.filter_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {

  if (length(x) == 1L) {
    return(x[[1L]])
  }

  if (!(length(unique(msaply(x, length))) == 1L &&
        length(unique(msaply(x, max))) == 1L &&
        length(unique(msaply(x, min))) == 1L)) {
    stop("Spectra differ in 'w.length' vector.")
  }

  if (getMultipleWl(x[[1]]) > 1L) {
    stop("Multiple spectra in long form not supported.")
  }

  is.A <- unique(msaply(x, is_absorbance_based))
  if (length(is.A) > 1L) {
    warning("Some spectra contain A and other Tfr: converting all to Tfr.")
    is.A <- FALSE
  }

  x <- A2T(x, action = "replace")

  num.spct <- length(x)
  w.length <- x[[1]][["w.length"]]
  M <- matrix(numeric(length(w.length) * num.spct),
              ncol = num.spct)

  for (i in seq_len(num.spct)) {
    M[ , i] <- x[[i]][["Tfr"]]
  }

  z <- apply(M, MARGIN = 1, FUN = mean, trim = trim, na.rm = na.rm)

  zz <- filter_spct(w.length = w.length,
                    Tfr = z,
                    Tfr.type = unique(msaply(x, getTfrType)),
                    multiple.wl = 1L)

  if (is.A) {
    T2A(zz)
  } else {
    zz
  }

}

#' @rdname mean
#'
#' @export
#'
mean.irrad_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {

  if (length(x) == 1L) {
    return(x[[1L]])
  }

  if (!(length(unique(msaply(x, length))) == 1L &&
        length(unique(msaply(x, max))) == 1L &&
        length(unique(msaply(x, min))) == 1L)) {
    stop("Spectra differ in 'w.length' vector.")
  }

  if (getMultipleWl(x[[1]]) > 1L) {
    stop("Multiple spectra in long form not supported.")
  }

  is.photon <- unique(msaply(x, is_photon_based))
  if (length(is.photon) > 1L) {
    warning("Some spectra contain s.q.irrad and others s.e.irrad: converting all to s.e.irrad.")
    is.photon <- FALSE
  }
  x <- q2e(x, action = "replace")

  num.spct <- length(x)
  w.length <- x[[1]][["w.length"]]
  M <- matrix(numeric(length(w.length) * num.spct),
              ncol = num.spct)

  for (i in seq_len(num.spct)) {
    M[ , i] <- x[[i]][["s.e.irrad"]]
  }

  z <- apply(M, MARGIN = 1, FUN = mean, trim = trim, na.rm = na.rm)

  zz <- irrad_spct(w.length = w.length,
                   s.e.irrad = z,
                   time.unit = unique(msaply(x, getTimeUnit)),
                   multiple.wl = 1L)

  if (is.photon) {
    e2q(zz)
  } else {
    zz
  }

}

#' @rdname mean
#'
#' @export
#'
mean.reflectance_mspct <- function(x, trim = 0, na.rm = FALSE, ...) {

  if (length(x) == 1L) {
    return(x[[1L]])
  }

  if (!(length(unique(msaply(x, length))) == 1L &&
        length(unique(msaply(x, max))) == 1L &&
        length(unique(msaply(x, min))) == 1L)) {
    stop("Spectra differ in 'w.length' vector.")
  }

  if (getMultipleWl(x[[1]]) > 1L) {
    stop("Multiple spectra in long form not supported.")
  }

  x <- q2e(x, action = "replace")

  num.spct <- length(x)
  w.length <- x[[1]][["w.length"]]
  M <- matrix(numeric(length(w.length) * num.spct),
              ncol = num.spct)

  for (i in seq_len(num.spct)) {
    M[ , i] <- x[[i]][["Rfr"]]
  }

  z <- apply(M, MARGIN = 1, FUN = mean, trim = trim, na.rm = na.rm)

  zz <- reflector_spct(w.length = w.length,
                       Rfr = z,
                       Rfr.type = unique(msaply(x, getRfrType)),
                       multiple.wl = 1L)

  zz

}
