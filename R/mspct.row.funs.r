#' rowwise functions for collections of spectra
#'
#' Private methods for collections of spectra. Applies a function at each wavelength
#' across all the spectra in the collection.
#'
#' @param x An R object. Currently this package defines methods for collections of
#'    spectral objects.
#' @param .fun An R function or a vector or list of functions.
#' @param col.name.tag character Vector of extensions to paste to default column
#'    name for the output from each of the functions. If col.name.tag[1] != "",
#'    this forces the return of an object of class \code{"generic_spct"}.
#' @param .fun.name character string used to set what.measured attribute.
#' @param ...	Arguments passed to .fun.
#'
#' @note Omission of NAs is done separately at each wavelength. Interpolation is
#'   not applied, so all spectra in \code{x} must share the same set of
#'   wavelengths.
#'
#' @seealso \code{\link[base]{mean}}
#'
#' @keywords internal
#'
#' @name rowwise
#'
rowwise_filter <- function(x, .fun, col.name.tag = "", .fun.name = "Summary of", ...) {

  # validate input
  stopifnot(length(.fun) == length(col.name.tag))

  if (!(length(unique(msaply(x, length))) == 1L &&
        length(unique(msaply(x, max))) == 1L &&
        length(unique(msaply(x, min))) == 1L)) {
    stop("Spectra differ in 'w.length' vector.")
  }

  if (getMultipleWl(x[[1]]) > 1L) {
    stop("Multiple spectra in long form not supported.")
  }

  # collect and check consistency of attributes
  is.Tfr <- unique(msaply(x, is_transmittance_based))
  if (length(is.Tfr) > 1L) {
    warning("Some spectra contain A and other Tfr: converting all to Tfr.")
    x <- A2T(x)
    is.Tfr <- TRUE
  }

  Tfr.type <- unique(msaply(x, getTfrType))
  if (length(Tfr.type) > 1L) {
    stop("Spectra differ in their Tfr.type: ", Tfr.type, ".")
  }

  # infer column name to use as input
  if (is.Tfr) {
    col.name <- "Tfr"
  } else {
    col.name <- "A"
  }

  # allocate memory
  num.spct <- length(x)
  w.length <- x[[1]][["w.length"]]
  M <- matrix(numeric(length(w.length) * num.spct),
              ncol = num.spct)

  # collect spectra into a matrix
  for (i in seq_len(num.spct)) {
    M[ , i] <- x[[i]][[col.name]]
  }

  # apply the functions and collect results
  z <- list(w.length)
  for (f in .fun) {
      z <- c(z, list(apply(M, MARGIN = 1, FUN = f, ...)))
  }
  col.name.out <- paste(col.name, col.name.tag, sep = "")
  names(z) <- c("w.length", col.name.out)
  zz <- tibble::as.tibble(z)

  # set class and attributes of spectrum to be returned
  if (col.name.tag[1] == "") {
    zz <- setFilterSpct(zz,
                        Tfr.type = Tfr.type,
                        multiple.wl = 1L)
  } else {
    zz <- setGenericSpct(zz,
                         multiple.wl = 1L)
  }
  setWhatMeasured(zz, paste(.fun.name, length(x), class(x[[1]])[1], "objects."))

  zz
}

#' @rdname rowwise
#'
#' @export
#'
rowwise_source <- function(x, .fun, col.name.tag = "", .fun.name = "Summary of", ...) {

  # validate input
  stopifnot(length(.fun) == length(col.name.tag))

  if (!(length(unique(msaply(x, length))) == 1L &&
        length(unique(msaply(x, max))) == 1L &&
        length(unique(msaply(x, min))) == 1L)) {
    stop("Spectra differ in 'w.length' vector.")
  }

  if (getMultipleWl(x[[1]]) > 1L) {
    stop("Multiple spectra in long form not supported.")
  }

  # collect and check consistency of attributes
  is.energy <- unique(msaply(x, is_energy_based))
  if (length(is.energy) > 1L) {
    warning("Some spectra contain 's.q.irrad' and others 's.e.irrad': converting all to s.e.irrad.")
    x <- q2e(x)
    is.energy <- TRUE
  }

  time.unit <- unique(msaply(x, getTimeUnit))
  if (length(time.unit) > 1L) {
    stop("Spectra differ in their 'time.unit' attribute: ", time.unit, ".")
  }

  bswf.used <- unique(msaply(x, getBSWFUsed))
  if (length(bswf.used) > 1L) {
    stop("Spectra differ in their 'bswf.used' attribute: ", bswf.used, ".")
  }

  # infer column name to use as input
  if (is.energy) {
    col.name <- "s.e.irrad"
  } else {
    col.name <- "s.q.irrad"
  }

  # allocate memory
  num.spct <- length(x)
  w.length <- x[[1]][["w.length"]]
  M <- matrix(numeric(length(w.length) * num.spct),
              ncol = num.spct)

  # collect spectra into a matrix
  for (i in seq_len(num.spct)) {
    M[ , i] <- x[[i]][[col.name]]
  }

  # apply the functions and collect results
  z <- list(w.length)
  for (f in .fun) {
    z <- c(z, list(apply(M, MARGIN = 1, FUN = f, ...)))
  }
  col.name.out <- paste(col.name, col.name.tag, sep = "")
  names(z) <- c("w.length", col.name.out)
  zz <- tibble::as.tibble(z)

  # set class and attributes of spectrum to be returned
  if (col.name.tag[1] == "") {
    zz <- setSourceSpct(zz,
                        time.unit = time.unit,
                        bswf.used = bswf.used,
                        multiple.wl = 1L)
  } else {
    zz <- setGenericSpct(zz,
                         multiple.wl = 1L)
  }
  setWhatMeasured(zz, paste(.fun.name, length(x), class(x[[1]])[1], "objects."))

  zz
}

#' @rdname rowwise
#'
#' @export
#'
rowwise_response <- function(x, .fun, col.name.tag = "", .fun.name = "Summary of", ...) {

  # validate input
  stopifnot(length(.fun) == length(col.name.tag))

  if (!(length(unique(msaply(x, length))) == 1L &&
        length(unique(msaply(x, max))) == 1L &&
        length(unique(msaply(x, min))) == 1L)) {
    stop("Spectra differ in 'w.length' vector.")
  }

  if (getMultipleWl(x[[1]]) > 1L) {
    stop("Multiple spectra in long form not supported.")
  }

  # collect and check consistency of attributes
  is.energy <- unique(msaply(x, is_energy_based))
  if (length(is.energy) > 1L) {
    warning("Some spectra contain 's.q.irrad' and others 's.e.irrad': converting all to s.e.irrad.")
    x <- q2e(x)
    is.energy <- TRUE
  }

  time.unit <- unique(msaply(x, getTimeUnit))
  if (length(time.unit) > 1L) {
    stop("Spectra differ in their 'time.unit' attribute: ", time.unit, ".")
  }

  # infer column name to use as input
  if (is.energy) {
    col.name <- "s.e.response"
  } else {
    col.name <- "s.q.response"
  }

  # allocate memory
  num.spct <- length(x)
  w.length <- x[[1]][["w.length"]]
  M <- matrix(numeric(length(w.length) * num.spct),
              ncol = num.spct)

  # collect spectra into a matrix
  for (i in seq_len(num.spct)) {
    M[ , i] <- x[[i]][[col.name]]
  }

  # apply the functions and collect results
  z <- list(w.length)
  for (f in .fun) {
    z <- c(z, list(apply(M, MARGIN = 1, FUN = f, ...)))
  }
  col.name.out <- paste(col.name, col.name.tag, sep = "")
  names(z) <- c("w.length", col.name.out)
  zz <- tibble::as.tibble(z)

  # set class and attributes of spectrum to be returned
  if (col.name.tag[1] == "") {
    zz <- setResponseSpct(zz,
                          time.unit = time.unit,
                          multiple.wl = 1L)
  } else {
    zz <- setGenericSpct(zz,
                         multiple.wl = 1L)
  }
  setWhatMeasured(zz, paste(.fun.name, length(x), class(x[[1]])[1], "objects."))

  zz
}

#' @rdname rowwise
#'
rowwise_reflector <- function(x, .fun, col.name.tag = "", .fun.name = "Summary of", ...) {

  # validate input
  stopifnot(length(.fun) == length(col.name.tag))

  if (!(length(unique(msaply(x, length))) == 1L &&
        length(unique(msaply(x, max))) == 1L &&
        length(unique(msaply(x, min))) == 1L)) {
    stop("Spectra differ in 'w.length' vector.")
  }

  if (getMultipleWl(x[[1]]) > 1L) {
    stop("Multiple spectra in long form not supported.")
  }

  # collect and check consistency of attributes
  Rfr.type <- unique(msaply(x, getRfrType))
  if (length(Rfr.type) > 1L) {
    stop("Spectra differ in their Rfr.type: ", Rfr.type, ".")
  }

  # infer column name to use as input
  col.name <- "Rfr"

  # allocate memory
  num.spct <- length(x)
  w.length <- x[[1]][["w.length"]]
  M <- matrix(numeric(length(w.length) * num.spct),
              ncol = num.spct)

  # collect spectra into a matrix
  for (i in seq_len(num.spct)) {
    M[ , i] <- x[[i]][[col.name]]
  }

  # apply the functions and collect results
  z <- list(w.length)
  for (f in .fun) {
    z <- c(z, list(apply(M, MARGIN = 1, FUN = f, ...)))
  }
  col.name.out <- paste(col.name, col.name.tag, sep = "")
  names(z) <- c("w.length", col.name.out)
  zz <- tibble::as.tibble(z)

  # set class and attributes of spectrum to be returned
  if (col.name.tag[1] == "") {
    zz <- setReflectorSpct(zz,
                        Rfr.type = Rfr.type,
                        multiple.wl = 1L)
  } else {
    zz <- setGenericSpct(zz,
                         multiple.wl = 1L)
  }
  setWhatMeasured(zz, paste(.fun.name, length(x), class(x[[1]])[1], "objects."))

  zz
}
