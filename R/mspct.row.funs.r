#' rowwise functions for collections of spectra
#'
#' Private methods for collections of spectra. Applies a function at each wavelength
#' across all the spectra in the collection.
#'
#' @param x An R object. Currently this package defines methods for collections of
#'    spectral objects.
#' @param .fun An R function or a list of functions.
#' @param col.name.tag character Vector of extensions to paste to default column
#'    name for the output from each of the functions. If \code{col.name.tag[1] != ""},
#'    this forces the return of an object of class \code{"generic_spct"}.
#' @param .fun.name character string used to set what.measured attribute.
#' @param w.length.out numeric vector of wavelengths (nanometres).
#' @param ...	Arguments passed to .fun.
#'
#' @note Omission of NAs is done separately at each wavelength. Interpolation is
#'   applied, so spectra in \code{x} do not need to share the same set of
#'   wavelengths. When defining
#'   new public functions using these utility functions make sure to return data
#'   that is valid for the class of spectral object returned!!
#'
#'   Objects of classes \code{raw_spct} and \code{cps_spct} can contain data
#'   from multiple scans acquired with the same instrument, thus at the same set
#'   of wave lengths.
#'
#' @keywords internal
#'
#' @name rowwise
#'
rowwise_filter <-
  function(x,
           .fun,
           col.name.tag = "",
           .fun.name = "Summary of",
           w.length.out = NULL,
           ...) {
    # we accept both function objects and lists of functions as input, but we
    # convert function arguments to lists of length 1..
    if (is.function(.fun)) {
      .fun <- list(.fun)
    }

    # validate input
    stopifnot(length(.fun) == length(col.name.tag))

    # ensure rowwise computations are possible
    x <- make_wl_consistent(x = x, w.length.out = w.length.out)

    # collect and check consistency of quantities
    is.Tfr <- all(msaply(x, is_transmittance_based))
    is.Afr <- all(msaply(x, is_absorptance_based))
    is.A <- all(msaply(x, is_absorbance_based))
    is.mixed <- !(is.Tfr || is.Afr || is.A)
    if (is.mixed) {
      warning("Not all spectra contain the same quantity: converting all to Tfr.")
      x <- any2T(x)
      is.Tfr <- TRUE
    }

    # infer column name to use as input
    if (is.Tfr) {
      col.name <- "Tfr"
    } else if (is.A) {
      col.name <- "A"
    } else {
      col.name <- "Afr"
    }

    Tfr.type <- unique(msaply(x, getTfrType))
    if (length(Tfr.type) > 1L) {
      stop("Spectra differ in their Tfr.type: ", Tfr.type, ".")
    }

    # allocate memory
    num.spct <- length(x)
    w.length <- x[[1]][["w.length"]]
    M <- matrix(numeric(length(w.length) * num.spct),
                ncol = num.spct)

    # collect spectra into a matrix
    for (i in seq_len(num.spct)) {
      M[, i] <- x[[i]][[col.name]]
    }

    # apply the functions and collect results
    z <- list(w.length)
    for (f in .fun) {
      z <- c(z, list(apply(
        M, MARGIN = 1, FUN = f, ...
      )))
    }
    col.name.out <- paste(col.name, col.name.tag, sep = "")
    names(z) <- c("w.length", col.name.out)
    zz <- tibble::as_tibble(z)

    # set class and attributes of spectrum to be returned
    if (col.name.tag[1] == "") {
      zz <- setFilterSpct(zz,
                          Tfr.type = Tfr.type,
                          multiple.wl = 1L)
    } else {
      zz <- setGenericSpct(zz,
                           multiple.wl = 1L)
    }
    what.measured <- what_measured(x, simplify = TRUE)
    what_measured(zz) <- paste(.fun.name, length(x), class(x[[1]])[1], "objects.\n",
                              what.measured)
    how_measured(zz) <- how_measured(x, simplify = TRUE)
    when_measured(zz) <- when_measured(x, simplify = TRUE)
    where_measured(zz) <- where_measured(x, simplify = TRUE)

    zz
  }

#' @rdname rowwise
#'
#' @export
#'
rowwise_source <-
  function(x,
           .fun,
           col.name.tag = "",
           .fun.name = "Summary of",
           w.length.out = NULL,
           ...) {
    # we accept both function objects and lists of functions as input, but we
    # convert function arguments to lists of length 1..
    if (is.function(.fun)) {
      .fun <- list(.fun)
    }

    # validate input
    stopifnot(length(.fun) == length(col.name.tag))

    # ensure rowise computtaions are possible
    x <- make_wl_consistent(x = x, w.length.out = w.length.out)

    # collect and check consistency of attributes
    is.energy <- unique(msaply(x, is_energy_based))
    if (length(is.energy) > 1L) {
      warning(
        "Some spectra contain 's.q.irrad' and others 's.e.irrad': converting all to s.e.irrad."
      )
      x <- q2e(x)
      is.energy <- TRUE
    }

    time.unit <- unique(msaply(x, getTimeUnit))
    if (length(time.unit) > 1L) {
      stop("Spectra differ in their 'time.unit' attribute: ",
           time.unit,
           ".")
    }

    bswf.used <- unique(msaply(x, getBSWFUsed))
    if (length(bswf.used) > 1L) {
      stop("Spectra differ in their 'bswf.used' attribute: ",
           bswf.used,
           ".")
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
      M[, i] <- x[[i]][[col.name]]
    }

    # apply the functions and collect results
    z <- list(w.length)
    for (f in .fun) {
      z <- c(z, list(apply(
        M, MARGIN = 1, FUN = f, ...
      )))
    }
    col.name.out <- paste(col.name, col.name.tag, sep = "")
    names(z) <- c("w.length", col.name.out)
    zz <- tibble::as_tibble(z)

    # set class and attributes of spectrum to be returned
    if (col.name.tag[1] == "") {
      zz <- setSourceSpct(
        zz,
        time.unit = time.unit,
        bswf.used = bswf.used,
        multiple.wl = 1L
      )
    } else {
      zz <- setGenericSpct(zz,
                           multiple.wl = 1L)
    }
    what.measured <- what_measured(x, simplify = TRUE)
    what_measured(zz) <- paste(.fun.name, length(x), class(x[[1]])[1], "objects.\n",
                               what.measured)
    how_measured(zz) <- how_measured(x, simplify = TRUE)
    when_measured(zz) <- when_measured(x, simplify = TRUE)
    where_measured(zz) <- where_measured(x, simplify = TRUE)

    zz
  }

#' @rdname rowwise
#'
#' @export
#'
rowwise_response <-
  function(x,
           .fun,
           col.name.tag = "",
           .fun.name = "Summary of",
           w.length.out = NULL,
           ...) {
    # we accept both function objects and lists of functions as input, but we
    # convert function arguments to lists of length 1..
    if (is.function(.fun)) {
      .fun <- list(.fun)
    }

    # validate input
    stopifnot(length(.fun) == length(col.name.tag))

    # ensure rowise computtaions are possible
    x <- make_wl_consistent(x = x, w.length.out = w.length.out)

    # collect and check consistency of attributes
    is.energy <- unique(msaply(x, is_energy_based))
    if (length(is.energy) > 1L) {
      warning(
        "Some spectra contain 's.q.irrad' and others 's.e.irrad': converting all to s.e.irrad."
      )
      x <- q2e(x)
      is.energy <- TRUE
    }

    time.unit <- unique(msaply(x, getTimeUnit))
    if (length(time.unit) > 1L) {
      stop("Spectra differ in their 'time.unit' attribute: ",
           time.unit,
           ".")
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
      M[, i] <- x[[i]][[col.name]]
    }

    # apply the functions and collect results
    z <- list(w.length)
    for (f in .fun) {
      z <- c(z, list(apply(
        M, MARGIN = 1, FUN = f, ...
      )))
    }
    col.name.out <- paste(col.name, col.name.tag, sep = "")
    names(z) <- c("w.length", col.name.out)
    zz <- tibble::as_tibble(z)

    # set class and attributes of spectrum to be returned
    if (col.name.tag[1] == "") {
      zz <- setResponseSpct(zz,
                            time.unit = time.unit,
                            multiple.wl = 1L)
    } else {
      zz <- setGenericSpct(zz,
                           multiple.wl = 1L)
    }
    what.measured <- what_measured(x, simplify = TRUE)
    what_measured(zz) <- paste(.fun.name, length(x), class(x[[1]])[1], "objects.\n",
                               what.measured)
    how_measured(zz) <- how_measured(x, simplify = TRUE)
    when_measured(zz) <- when_measured(x, simplify = TRUE)
    where_measured(zz) <- where_measured(x, simplify = TRUE)

    zz
  }

#' @rdname rowwise
#'
rowwise_reflector <-
  function(x,
           .fun,
           col.name.tag = "",
           .fun.name = "Summary of",
           w.length.out = NULL,
           ...) {
    # we accept both function objects and lists of functions as input, but we
    # convert function arguments to lists of length 1..
    if (is.function(.fun)) {
      .fun <- list(.fun)
    }

    # validate input
    stopifnot(length(.fun) == length(col.name.tag))

    # ensure rowise computtaions are possible
    x <- make_wl_consistent(x = x, w.length.out = w.length.out)

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
      M[, i] <- x[[i]][[col.name]]
    }

    # apply the functions and collect results
    z <- list(w.length)
    for (f in .fun) {
      z <- c(z, list(apply(
        M, MARGIN = 1, FUN = f, ...
      )))
    }
    col.name.out <- paste(col.name, col.name.tag, sep = "")
    names(z) <- c("w.length", col.name.out)
    zz <- tibble::as_tibble(z)

    # set class and attributes of spectrum to be returned
    if (col.name.tag[1] == "") {
      zz <- setReflectorSpct(zz,
                             Rfr.type = Rfr.type,
                             multiple.wl = 1L)
    } else {
      zz <- setGenericSpct(zz,
                           multiple.wl = 1L)
    }
    what.measured <- what_measured(x, simplify = TRUE)
    what_measured(zz) <- paste(.fun.name, length(x), class(x[[1]])[1], "objects.\n",
                               what.measured)
    how_measured(zz) <- how_measured(x, simplify = TRUE)
    when_measured(zz) <- when_measured(x, simplify = TRUE)
    where_measured(zz) <- where_measured(x, simplify = TRUE)

    zz
  }

#' @rdname rowwise
#'
rowwise_calibration <-
  function(x,
           .fun,
           col.name.tag = "",
           .fun.name = "Summary of",
           w.length.out = NULL,
           ...) {
    # we accept both function objects and lists of functions as input, but we
    # convert function arguments to lists of length 1..
    if (is.function(.fun)) {
      .fun <- list(.fun)
    }

    # validate input
    stopifnot(length(.fun) == length(col.name.tag))

    # ensure rowise computtaions are possible
    x <- make_wl_consistent(x = x, w.length.out = w.length.out)

    # infer column name to use as input
    col.name <- "irrad.mult"

    # allocate memory
    num.spct <- length(x)
    w.length <- x[[1]][["w.length"]]
    M <- matrix(numeric(length(w.length) * num.spct),
                ncol = num.spct)

    # collect spectra into a matrix
    for (i in seq_len(num.spct)) {
      M[, i] <- x[[i]][[col.name]]
    }

    # apply the functions and collect results
    z <- list(w.length)
    for (f in .fun) {
      z <- c(z, list(apply(
        M, MARGIN = 1, FUN = f, ...
      )))
    }
    col.name.out <- paste(col.name, col.name.tag, sep = "")
    names(z) <- c("w.length", col.name.out)
    zz <- tibble::as_tibble(z)

    # set class and attributes of spectrum to be returned
    zz <- setCalibrationSpct(zz, multiple.wl = 1L)

    what.measured <- what_measured(x, simplify = TRUE)
    what_measured(zz) <- paste(.fun.name, length(x), class(x[[1]])[1], "objects.\n",
                               what.measured)
    how_measured(zz) <- how_measured(x, simplify = TRUE)
    when_measured(zz) <- when_measured(x, simplify = TRUE)
    where_measured(zz) <- where_measured(x, simplify = TRUE)

    zz
  }

#' @rdname rowwise
#'
rowwise_cps <-
  function(x,
           .fun,
           col.name.tag = "",
           .fun.name = "Summary of",
           w.length.out = NULL,
           ...) {
    # we accept both function objects and lists of functions as input, but we
    # convert function arguments to lists of length 1..
    if (is.function(.fun)) {
      .fun <- list(.fun)
    }

    # validate input
    stopifnot(length(.fun) == length(col.name.tag))

    # ensure rowise computtaions are possible
    x <- make_wl_consistent(x = x, w.length.out = w.length.out)

    # infer column name to use as input
    col.name = "cps"
    have.cps.colum <- sapply(x,
                            FUN = function(y) {
                              sum(grepl("^cps$", colnames(y))) == 1L
                            })
    if (!all(have.cps.colum)) {
      message("Columns differ among spectra.")
      zz <- cps_spct(w.length = x[[1]][["w.length"]], cps = NA_real_)
    } else {
      # allocate memory
      num.spct <- length(x)
      w.length <- x[[1]][["w.length"]]
      M <- matrix(numeric(length(w.length) * num.spct),
                  ncol = num.spct)

      # collect spectra into a matrix
      for (i in seq_len(num.spct)) {
        M[, i] <- x[[i]][[col.name]]
      }

      # apply the functions and collect results
      z <- list(w.length)
      for (f in .fun) {
        z <- c(z, list(apply(
          M, MARGIN = 1, FUN = f, ...
        )))
      }
      col.name.out <- paste(col.name, col.name.tag, sep = "")
      names(z) <- c("w.length", col.name.out)
      zz <- tibble::as_tibble(z)

      # set class and attributes of spectrum to be returned
      zz <- setCpsSpct(zz, multiple.wl = 1L)
    }

    what.measured <- what_measured(x, simplify = TRUE)
    what_measured(zz) <- paste(.fun.name, length(x), class(x[[1]])[1], "objects.\n",
                               what.measured)
    how_measured(zz) <- how_measured(x, simplify = TRUE)
    when_measured(zz) <- when_measured(x, simplify = TRUE)
    where_measured(zz) <- where_measured(x, simplify = TRUE)

    zz
  }

#' @rdname rowwise
#'
rowwise_raw <-
  function(x,
           .fun,
           col.name.tag = "",
           .fun.name = "Summary of",
           w.length.out = NULL,
           ...) {
    # we accept both function objects and lists of functions as input, but we
    # convert function arguments to lists of length 1..
    if (is.function(.fun)) {
      .fun <- list(.fun)
    }

    # validate input
    stopifnot(length(.fun) == length(col.name.tag))

    # ensure rowise computtaions are possible
    x <- make_wl_consistent(x = x, w.length.out = w.length.out)

    # infer column name to use as input
    col.name = "counts"
    have.counts.colum <- sapply(x,
                                FUN = function(y) {
                                  sum(grepl("^counts$", colnames(y))) == 1L
                                })
    if (!all(have.counts.colum)) {
      message("Columns differ among spectra.")
      zz <- cps_spct(w.length = x[[1]][["w.length"]], cps = NA_real_)
    } else {
      # allocate memory
      num.spct <- length(x)
      w.length <- x[[1]][["w.length"]]
      M <- matrix(numeric(length(w.length) * num.spct),
                  ncol = num.spct)

      # collect spectra into a matrix
      for (i in seq_len(num.spct)) {
        M[, i] <- x[[i]][[col.name]]
      }

      # apply the functions and collect results
      z <- list(w.length)
      for (f in .fun) {
        z <- c(z, list(apply(
          M, MARGIN = 1, FUN = f, ...
        )))
      }
      col.name.out <- paste(col.name, col.name.tag, sep = "")
      names(z) <- c("w.length", col.name.out)
      zz <- tibble::as_tibble(z)

      # set class and attributes of spectrum to be returned
      zz <- setRawSpct(zz, multiple.wl = 1L)
    }

    what.measured <- what_measured(x, simplify = TRUE)
    what_measured(zz) <- paste(.fun.name, length(x), class(x[[1]])[1], "objects.\n",
                               what.measured)
    how_measured(zz) <- how_measured(x, simplify = TRUE)
    when_measured(zz) <- when_measured(x, simplify = TRUE)
    where_measured(zz) <- where_measured(x, simplify = TRUE)

    zz
  }

#' Make wavelengths consistent
#'
#' Ensure wavelengths in all spectra in a \code{mspct} object are consistent,
#' filling head and tail of spectral quantities when needed with NAs, and using
#' interpolation if needed. If individual \code{spct} objects contain multiple
#' spectra in long form these are split into multiple \code{spct} objects, each
#' containing data for one spectrum.
#'
#' @param x a generic_spct object
#' @param w.length.out numeric vector of wavelengths (nanometres)
#'
#' @keywords internal
#'
make_wl_consistent <- function(x,
                               w.length.out = NULL) {
  # split multiple spectra in long form
  if (any(msaply(x, .fun = getMultipleWl) > 1L)) {
    x <- subset2mspct(x)
    message("Multiple spectra in long form have been split.")
  }

  if (length(x) <= 1L) {
    # nothing else to do
    return(x)
  }

  # interpolate wavelengths to make them consistent
  w.lengths.differ <-
    !(length(unique(msaply(x, .fun = nrow))) == 1L &&
        length(unique(msaply(x, .fun = wl_max))) == 1L &&
        length(unique(msaply(x, .fun = wl_min))) == 1L)

  if (w.lengths.differ || !is.null(w.length.out)) {
    if (is.null(w.length.out)) {
      if (length(x) == 2) {
        # we use the union of wavelengths as binary operators do
        w.length.out <-
          sort(union(x[[1]][["w.length"]], x[[1]][["w.length"]]))
      } else {
        # we avoid bloating of the data while protection the wl resolution
        wl.ranges <- range(x, na.rm = TRUE)
        wl.range <- numeric(2)
        wl.range[1] <- min(wl.ranges["min.wl"])
        wl.range[2] <- max(wl.ranges["max.wl"])
        max.length <- max(msaply(x, .fun = length))
        stepsize <- min(msaply(x, .fun = function(x) { stats::median(diff(x[["w.length"]])) }))
        length.out <- max(max.length, diff(wl.range) / stepsize)
        w.length.out <-
          seq.int(from = wl.range[1], to = wl.range[2], length.out = length.out)
      }
    }
    x <- interpolate_mspct(x, w.length.out = w.length.out, fill = NA)
  }
  x
}
