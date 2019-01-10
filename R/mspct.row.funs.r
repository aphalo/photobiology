#' rowwise functions for collections of spectra
#'
#' Private methods for collections of spectra. Applies a function at each wavelength
#' across all the spectra in the collection.
#'
#' @param x An R object. Currently this package defines methods for collections of
#'    spectral objects.
#' @param .fun An R function.
#' @param ...	Arguments passed to .fun.
#'
#' @note Omission of NAs is done separately at each wavelength. Interpolation is
#'   not applied, so all spectra in \code{x} must share the same set of
#'   wavelengths.
#'
#' @seealso \code{\link[base]{mean}}
#'
#' @name rowwise
#'

rowwise_filter <- function(x, .fun, ...) {

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

  z <- apply(M, MARGIN = 1, FUN = .fun, ...)

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

#' @rdname rowwise
#'
#' @export
#'
rowwise_irrad <- function(x, .fun, ...) {

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

  z <- apply(M, MARGIN = 1, FUN = .fun, ...)

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

#' @rdname rowwise
#'
rowwise_Rfr <- function(x, .fun, ...) {

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

  z <- apply(M, MARGIN = 1, FUN = .fun, ...)

  zz <- reflector_spct(w.length = w.length,
                       Rfr = z,
                       Rfr.type = unique(msaply(x, getRfrType)),
                       multiple.wl = 1L)

  zz

}
