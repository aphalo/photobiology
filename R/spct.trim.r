#' Trim (or expand) head and/or tail of a spectrum
#'
#' Trim head and tail of a spectrum based on wavelength limits, interpolating
#' the values at the boundaries of the range. Trimming is needed for example to
#' remove short wavelength noise when the measured spectrum extends beyond the
#' known emission spectrum of the measured light source. Occasionally one may
#' want also to expand the wavelength range.
#'
#' @param spct an object of class "generic_spct".
#' @param range a numeric vector of length two, or any other object for which
#'   method range() will return a numeric vector of length two.
#' @param low.limit shortest wavelength to be kept (defaults to shortest
#'   w.length value).
#' @param high.limit longest wavelength to be kept (defaults to longest w.length
#'   value).
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param fill if fill==NULL then tails are deleted, otherwise tails or s.irrad
#'   are filled with the value of fill.
#' @param byref logical indicating if new object will be created by reference or
#'   by copy of spct.
#' @param verbose logical.
#'
#' @return a spectrum object or a collection of spectral objects of the same
#'   class as \code{x} with wavelength heads and tails clipped or extended.
#'
#' @note When expanding a spectrum, if fill==NULL, then expansion is not
#'   performed. Range can be "waveband" object, a numeric vector or a list of
#'   numeric vectors, or any other user-defined or built-in object for which
#'   \code{range()} returns a numeric vector of length two, that can be
#'   interpreted as wavelengths expressed in nm.
#' @family trim functions
#'
#' @export
#' @examples
#' trim_spct(sun.spct, low.limit=300)
#' trim_spct(sun.spct, low.limit=300, fill=NULL)
#' trim_spct(sun.spct, low.limit=300, fill=NA)
#' trim_spct(sun.spct, low.limit=300, fill=0.0)
#' trim_spct(sun.spct, range = c(300, 400))
#' trim_spct(sun.spct, range = c(300, NA))
#' trim_spct(sun.spct, range = c(NA, 400))
#'
trim_spct <- function(spct,
                      range = NULL,
                      low.limit = NULL, high.limit = NULL,
                      use.hinges = TRUE,
                      fill = NULL,
                      byref = FALSE,
                      verbose = getOption("photobiology.verbose") )
{
  if (nrow(spct) == 0) {
    return(spct)
  }
  stopifnot(is.generic_spct(spct))
  x <- spct
  num.spectra <- getMultipleWl(x)
  if (num.spectra != 1) {
    warning("Skipping trim operation as object contains ",
            num.spectra, " spectra")
    return(x)
  }
  if (is.null(use.hinges)) {
    use.hinges <- auto_hinges(spct[["w.length"]])
  }
  if (byref) {
    name <- substitute(spct)
  }
  # class_spct <- class(spct)
  if (is.null(range)) {
    range[1] <- ifelse(is.null(low.limit), NA, low.limit)
    range[2] <- ifelse(is.null(high.limit), NA, high.limit)
  }
  range <- normalize_range_arg(range, range(spct), trim = FALSE) # trim = is.null(fill)?
  low.limit <- range[1]
  high.limit <- range[2]
  trim.low <- low.limit > min(spct, na.rm = TRUE)
  trim.high <- high.limit < max(spct, na.rm = TRUE)
  if (trim.low && trim.high && high.limit - low.limit < 1e-7) {
    warning("When trimming, 'range' must be a finite wavelength interval > 1E-7 nm")
    return(spct[FALSE, ]) # returns a spct object with nrow equal to zero
  }
  names.spct <- names(spct)
  names.data <- names.spct[names.spct != "w.length"]

  # check whether we should expand the low end
  low.end <- min(spct, na.rm = TRUE)
  if (!trim.low && low.end > low.limit) {
    if (!is.null(fill)) {
      # expand short tail
      low.tail.length <-  trunc(low.end - low.limit) + ifelse(use.hinges, 2, 1)
      low.tail.w.length <- seq(from = low.limit,
                               to = ifelse(use.hinges, low.end - 1e-12, low.end - 1),
                               length.out = low.tail.length)
      spct.top <- tibble::tibble(w.length = low.tail.w.length)
      for (data.col in names.data) {
        col.class <- class(spct[[data.col]])[1]
        if ("numeric" %in% col.class) {
          spct.top[[data.col]] <- fill
        } else if ("character" %in% col.class) {
          if (length(unique(spct.top[[data.col]])) == 1L) {
            spct.top[[data.col]] <- spct[1, data.col]
          } else {
            spct.top[[data.col]] <- NA_character_
          }
        } else if ("factor" %in% col.class) {
          if (length(levels(spct.top[[data.col]])) == 1L) {
            spct.top[[data.col]] <- spct[1, data.col]
          } else {
            spct.top[[data.col]] <- NA
          }
        } else {
          spct.top[[data.col]] <- NA
        }
      }
      spct <- plyr::rbind.fill(list(spct.top, spct))
      spct <- tibble::as_tibble(spct)
      setGenericSpct(spct)
      low.end <- min(spct)
    } else {
        if (verbose && (low.end - low.limit) > 0.01) {
          warning("Not trimming short end as low.limit is outside spectral data range.")
        }
      trim.low <- FALSE
    }
  }

  # check whether we should expand the high end
  high.end <- max(spct, na.rm = TRUE)
  if (!trim.high && high.end < high.limit) {
    if (!is.null(fill)) {
      # expand short tail
      high.tail.length <- trunc(high.limit - high.end) + ifelse(use.hinges, 2, 1)
      high.tail.w.length <- seq(from = ifelse(use.hinges, high.end + 1e-12, high.end + 1),
                                to = high.limit,
                                length.out = high.tail.length)
      spct.bottom <- tibble::tibble(w.length = high.tail.w.length)
      for (data.col in names.data) {
        col.class <- class(spct[[data.col]])[1]
        if ("numeric" %in% col.class) {
          spct.bottom[[data.col]] <- fill
        } else if ("character" %in% col.class) {
          if (length(unique(spct.bottom[[data.col]])) == 1L) {
            spct.bottom[[data.col]] <- spct[1, data.col]
          } else {
            spct.bottom[[data.col]] <- NA_character_
          }
        } else if ("factor" %in% col.class) {
          if (length(levels(spct.bottom[[data.col]])) == 1L) {
            spct.bottom[[data.col]] <- spct[1, data.col]
          } else {
            spct.bottom[[data.col]] <- NA
          }
        } else {
          spct.bottom[[data.col]] <- NA
        }
      }
      spct <- plyr::rbind.fill(list(spct, spct.bottom))
      spct <- tibble::as_tibble(spct)
      setGenericSpct(spct)
      low.end <- max(spct)
    } else {
      # give a warning only if difference is > 0.01 nm
      if (verbose && (high.limit - high.end) > 0.01) {
        warning("Not trimming long end as high.limit is outside spectral data range.")
      }
     trim.high <- FALSE
    }
  }

  # insert hinges
  if (use.hinges) {
    hinges <- NULL
    if (trim.low) {
      hinges <- c(hinges, low.limit - 1e-12, low.limit)
    }
    if (trim.high) {
      hinges <- c(hinges, high.limit - 1e-12, high.limit)
    }
    spct <- insert_spct_hinges(spct, hinges)
  }

  if (trim.low && trim.high) {
    within.selector <- with(spct, w.length >= low.limit & w.length < high.limit)
  } else if (trim.low) {
    within.selector <- with(spct, w.length >= low.limit)
  } else if (trim.high) {
    within.selector <- with(spct, w.length < high.limit)
  } else {
    within.selector <- TRUE
  }
  if (is.null(fill)) {
    spct <- spct[within.selector, ]
  } else {
    for (data.col in names.data) {
      if (is.numeric(spct[[data.col]])) {
       spct[!within.selector, data.col] <- fill
      }
    }
  }
  # we still need to copy all attributes as when we use row bind to
  # data.frames to extend the spectra.
  spct <- copy_attributes(x, spct, copy.class = TRUE)
  check_spct(spct)
  if (byref && is.name(name)) {
    name <- as.character(name)
    assign(name, spct, parent.frame(), inherits = TRUE)
  }
  spct
}


#' @rdname trim_spct
#'
#' @param mspct an object of class "generic_mspct"
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
trim_mspct <- function(mspct,
                       range = NULL,
                       low.limit = NULL,
                       high.limit = NULL,
                       use.hinges = TRUE,
                       fill = NULL,
                       byref = FALSE,
                       verbose = getOption("photobiology.verbose"),
                       .parallel = FALSE,
                       .paropts = NULL) {

  if (!length(mspct)) return(mspct) # class of x in no case changes

  name <- substitute(mspct)

  z <- msmsply(mspct = mspct,
               .fun = trim_spct,
               range = range,
               low.limit = low.limit,
               high.limit = high.limit,
               use.hinges = use.hinges,
               fill = fill,
               byref = FALSE,
               verbose = verbose,
               .parallel = .parallel,
               .paropts = .paropts)

  z <- copy_attributes(x = mspct, y = z, copy.class = TRUE) # msmsply strips attributes!!

  if (byref & is.name(name)) {
    name <- as.character(name)
    assign(name, z, parent.frame(), inherits = TRUE)
  }
  z
}

#' @rdname trim_spct
#'
#' @export
#'
trim2overlap <- function(mspct,
                         use.hinges = TRUE,
                         verbose = getOption("photobiology.verbose"),
                         .parallel = FALSE,
                         .paropts = NULL) {
  stopifnot(is.any_mspct(mspct))
  if (length(mspct) < 2) {
    # nothing to do
    return(mspct)
  }
  ranges <- msdply(mspct, range,
                   .parallel = .parallel,
                   .paropts = .paropts)
  range <- with(ranges, c(max(min.wl), min(max.wl)))
  trim_mspct(mspct,
             range = range,
             use.hinges = use.hinges,
             fill = NULL,
             verbose = verbose,
             .parallel = .parallel,
             .paropts = .paropts)
}

#' @rdname trim_spct
#'
#' @export
#'
extend2extremes <- function(mspct,
                            use.hinges = TRUE,
                            fill = NA,
                            verbose = getOption("photobiology.verbose"),
                            .parallel = FALSE,
                            .paropts = NULL) {
  stopifnot(is.any_mspct(mspct))
  if (length(mspct) < 2) {
    # nothing to do
    return(mspct)
  }
  ranges <- msdply(mspct,
                   .fun = range,
                   .parallel = .parallel,
                   .paropts = .paropts)
  range <- with(ranges, c(min(min.wl), max(max.wl)))
  trim_mspct(mspct,
             range = range,
             use.hinges = use.hinges,
             fill = fill,
             verbose = verbose,
             .parallel = .parallel,
             .paropts = .paropts)
}

#' Trim head and/or tail of a spectrum
#'
#' Trim head and tail of a spectrum based on wavelength limits, with
#' interpolation at range boundaries used by default. Expansion is also
#' possible.
#'
#' @param x an R object.
#' @param range a numeric vector of length two, or any other object for which
#'   function range() will return two.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param fill if \code{fill == NULL} then tails are deleted, otherwise tails
#'   are filled with the value of fill.
#' @param ... ignored (possibly used by derived methods).
#'
#' @return A copy of \code{x}, usually trimmed or expanded to a different
#'   length, either shorter or longer. Possibly with some of the original
#'   spectral data values replaced with \code{fill}.
#'
#' @note By default the \code{w.length} values for the first and last rows
#'   in the returned object are the values supplied as \code{range}.
#'
#' @family trim functions
#' @export
#' @examples
#' trim_wl(sun.spct, range = c(400, 500))
#' trim_wl(sun.spct, range = c(NA, 500))
#' trim_wl(sun.spct, range = c(400, NA))
#'
trim_wl <- function(x, range, use.hinges, fill, ...) UseMethod("trim_wl")

#' @describeIn trim_wl Default for generic function
#'
#' @export
#'
trim_wl.default <- function(x, range, use.hinges, fill, ...) {
  warning("'trim_wl' is not defined for objects of class ", class(x)[1])
  x
}

#' @describeIn trim_wl Trim an object of class "generic_spct" or derived.
#'
#' @export
#'
trim_wl.generic_spct <- function(x,
                                 range = NULL,
                                 use.hinges = TRUE,
                                 fill = NULL,
                                 ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(x) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    # call method on the collection
    z <- trim_wl(x = mspct,
                 range = range,
                 use.hinges = use.hinges,
                 fill = fill,
                 ...)
    return(rbindspct(z, idfactor = FALSE))
  }

  if (is.null(range)) {
    return(x)
  }
  trim_spct(spct = x,
            range = range,
            use.hinges = use.hinges,
            fill = fill,
            byref = FALSE,
            verbose = getOption("photobiology.verbose") )
}

#' @describeIn trim_wl  Trim an object of class "generic_mspct" or derived.
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
trim_wl.generic_mspct <- function(x,
                                  range = NULL,
                                  use.hinges = TRUE,
                                  fill = NULL, ...,
                                  .parallel = FALSE,
                                  .paropts = NULL) {

  if (!length(x)) return(x) # class of x in no case changes

  x <- subset2mspct(x) # expand long form spectra within collection

  if (is.null(range)) {
    return(x)
  }
  trim_mspct(mspct = x,
             range = range,
             use.hinges = use.hinges,
             fill = fill,
             byref = FALSE,
             verbose = getOption("photobiology.verbose"),
             .parallel = .parallel,
             .paropts = .paropts)
}

#' @describeIn trim_wl Trim an object of class "waveband".
#' @param trim logical (default is TRUE which trims the wavebands at the
#'   boundary, while FALSE discards wavebands that are partly off-boundary).
#'
#' @note trim_wl when applied to waveband objects always inserts hinges when
#'   trimming.
#'
#' @export
#'
trim_wl.waveband <- function(x,
                             range = NULL,
                             use.hinges = TRUE,
                             fill = NULL,
                             trim = getOption("photobiology.waveband.trim",
                                              default = TRUE),
                             ...) {
  if (is.null(range)) {
    return(x)
  }
  trim_waveband(w.band = x,
                range = range,
                trim = trim,
                use.hinges = use.hinges)
}

#' @describeIn trim_wl Trim a list (of "waveband" objects).
#'
#' @note trim_wl when applied to waveband objects always inserts hinges when
#'   trimming.
#'
#' @export
#'
trim_wl.list <- function(x,
                         range = NULL,
                         use.hinges = TRUE,
                         fill = NULL,
                         trim = getOption("photobiology.waveband.trim",
                                          default = TRUE),
                         ...) {

  if (!length(x)) return(x) # class of x in no case changes

  stopifnot(all(sapply(x, is.waveband)))
  if (is.null(range)) {
    return(x)
  }
  trim_waveband(w.band = x,
                range = range,
                trim = trim,
                use.hinges = use.hinges)
}

#' Clip head and/or tail of a spectrum
#'
#' Clip head and tail of a spectrum based on wavelength limits, no
#' interpolation used at range boundaries.
#'
#' @param x an R object.
#' @param range a numeric vector of length two, or any other object for which
#'   function \code{range()} will return range of wavelengths expressed in
#'   nanometres.
#' @param ... ignored (possibly used by derived methods).
#'
#' @return a spectrum object or a collection of spectral objects of the same
#'   class as \code{x} with wavelength heads and tails clipped.
#'
#' @note The condition tested is \code{wl >= range[1] & wl < (range[2] + 1e-13)}.
#'
#' @family trim functions
#' @export
#' @examples
#' clip_wl(sun.spct, range = c(400, 500))
#' clip_wl(sun.spct, range = c(NA, 500))
#' clip_wl(sun.spct, range = c(400, NA))
#'
clip_wl <- function(x, range, ...) UseMethod("clip_wl")

#' @describeIn clip_wl Default for generic function
#'
#' @export
#'
clip_wl.default <- function(x, range, ...) {
  warning("'clip_wl' is not defined for objects of class ", class(x)[1])
  x
}

#' @describeIn clip_wl Clip an object of class "generic_spct" or derived.
#'
#' @export
#'
clip_wl.generic_spct <- function(x, range = NULL, ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(x) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    # call method on the collection
    z <- trim_wl(x = mspct,
                 range = range,
                 ...)
    return(rbindspct(z, idfactor = FALSE))
  }

  if (is.null(range)) {
    return(x)
  }
  guard <- 1e-13
  stopifnot(is.generic_spct(x))
  stopifnot(!all(is.na(range)))
  if (is.numeric(range) && length(range) == 2) {
    if (is.na(range[1])) {
      x[x[["w.length"]] < range[2] + guard, ]
    } else if (is.na(range[2])) {
      x[x[["w.length"]] >= range[1], ]
    } else {
      x[x[["w.length"]] >= range[1] & x[["w.length"]] < range[2] + guard, ]
    }
  } else {
    range = range(range, na.rm = TRUE)
    x[x[["w.length"]] >= range[1] & x[["w.length"]] < range[2] + guard, ]
  }
}

#' @describeIn clip_wl  Clip an object of class "generic_mspct" or derived.
#'
#' @export
#'
clip_wl.generic_mspct <- function(x, range = NULL, ...) {
  if (!length(x)) return(x) # class of x in no case changes

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(mspct = x,
          .fun = clip_wl,
          range = range)
}

#' @describeIn clip_wl Clip an object of class "waveband".
#'
#' @export
#'
clip_wl.waveband <- function(x, range = NULL, ...) {
  if (is.null(range)) {
    return(x)
  }
  trim_waveband(w.band = x,
                range = range,
                trim = FALSE)
}

#' @describeIn clip_wl Clip a list (of objects of class "waveband").
#'
#' @export
#'
clip_wl.list <- function(x, range = NULL, ...) {

  if (!length(x)) return(x) # class of x in no case changes

  stopifnot(all(sapply(x, is.waveband)))
  if (is.null(range)) {
    return(x)
  }
  trim_waveband(w.band = x,
                range = range,
                trim = FALSE)
}
