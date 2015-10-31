
# Apply -------------------------------------------------------------------

#' Multi-spct transform methods
#'
#' Apply a function or operator to a collection of spectra.
#'
#' @param mspct an object of class generic_mspct or a derived class
#' @param .fun a function
#' @param ... other arguments passed to .fun
#'
#' @return a collection of spectra in the case of \code{msmsply}
#'
#' @export
#'
msmsply <- function(mspct, .fun, ...) {
  stopifnot(is.any_mspct(mspct))
  mspct.class <- class(mspct)
  byrow <- attr(mspct, "mspct.byrow", exact = TRUE)
  dim <- dim(mspct)
  ncol <- ncol(mspct)
  # llply returns a matrix for classes derived from list
  #
  rmDerivedMspct(mspct)
  y <- plyr::llply(mspct, .fun, ...)

  stopifnot(length(y) == length(mspct))

  if (length(y) > 1) {
    result.class <- shared_member_class(y)[1]
  } else {
    result.class <- mspct.class[1]
  }
  stopifnot(length(result.class) == 1)

  generic_mspct(l = y,
                class = result.class,
                byrow = byrow,
                ncol = ncol,
                dim = dim)
}

#' @rdname  msmsply
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   mspct, if \code{NULL}, the default, a column is added only if all members
#'   of \code{mscpt} are named.
#' @param col.names character Names to be used for data columns.
#'
#' @return a data frame in the case of \code{msdply}
#'
#' @export
#'
msdply <- function(mspct, .fun, ..., idx = NULL, col.names = NULL) {
  stopifnot(is.any_mspct(mspct))

  if ( (is.logical(idx) && idx) ||
       (is.null(idx) && !any(is.null(names(mspct)))) ) {
    .idx <- "spct.idx"
  } else if (is.logical(idx) && !idx) {
    .idx <- NULL
  } else {
    .idx <- idx
  }

  z <- plyr::ldply(.data = mspct,
                   .fun = .fun,
                   ...,
                   .id = .idx )

  f.name <- as.character( substitute(.fun))

  if (f.name %in% c("min",
                    "max",
                    "range",
                    "spread",
                    "midpoint",
                    "stepsize",
                    "getWhenMeasured",
                    "getWhereMeasured")) {
    qty.names <- switch(
      f.name,
      min = "min.wl",
      max = "max.wl",
      range = c("min.wl", "max.wl"),
      spread = "spread.wl",
      midpoint = "midpoint.wl",
      stepsize = c("min.step.wl", "max.step.wl"),
      getWhenMeasured = "when.measured",
      getWhereMeasured = NULL
    )
  } else if (!is.null(col.names) &&
             !any(col.names == "") &&
             !any(is.na(col.names)) &&
             length(col.names) == length(names(z)) - 1) {
    qty.names <- col.names
#  } else if (any(c("total", "mean", "contrib", "particip") %in% tolower(names(z)))) {
  } else {# make new names using function name
    if (is.null(.idx)) qty.names <- names(z) else qty.names <- names(z)[-1]
    qty.names <- paste(f.name,
                       gsub(" ", "", qty.names),
                       sep = "_")
#    qty.names <- NULL
  }

  if (!is.null(qty.names)) {
    if (is.null(.idx)) {
      names(z) <- qty.names
    } else {
      names(z)[-1] <- qty.names
    }
  }

  comment(z) <- paste("Applied function: '", f.name, "'.\n", sep = "", comment(mspct))

  mspct.nrow <- nrow(mspct)
  mspct.ncol <- ncol(mspct)
  mspct.byrow <- attr(mspct, "mspct.byrow", exact = TRUE)
  if (is.null(mspct.byrow)) {
    mspct.nrow <- FALSE
  }
  if (mspct.ncol > 1) {
    if (mspct.byrow) {
      z$col <- rep(1:mspct.ncol, mspct.nrow)
      z$row <- rep(1:mspct.nrow, rep(mspct.ncol, mspct.nrow))
    } else {
      z$col <- rep(1:mspct.ncol, rep(mspct.nrow, mspct.ncol))
      z$row <- rep(1:mspct.nrow, mspct.ncol)
    }
  }
  dplyr::as_data_frame(z)
}

#' @rdname  msmsply
#'
#' @return a list in the case of \code{mslply}
#'
#' @export
#'
mslply <- function(mspct, .fun, ...) {
  stopifnot(is.any_mspct(mspct))

  # llply returns a matrix for classes derived from list
  #
  rmDerivedMspct(mspct)
  z <- plyr::llply(.data = mspct,
                   .fun = .fun,
                   ...)

  names(z) <- names(mspct)

  f.name <- as.character( substitute(.fun))

  comment(z) <- paste("Applied function: '", f.name, "'.\n", sep = "", comment(mspct))

  z
}


#' @rdname  msmsply
#'
#' @param .drop should extra dimensions of length 1 in the output be dropped,
#'   simplifying the output. Defaults to TRUE
#' @return an array in the case of \code{msaply}
#'
#' @export
#'
msaply <- function(mspct, .fun, ..., .drop = TRUE) {
  stopifnot(is.any_mspct(mspct))

  # As many of our summary functions return nuneric values with names and other
  # attributes they need to be removed for dply::lapply to accept them.
  .ffun <- function(mspct, ...) {
    z <- .fun(mspct, ...)
    if (is.numeric(z)) {
      z <- as.numeric(z)
    }
    z
  }

  z <- plyr::laply(.data = mspct,
                   .fun = .ffun,
                   ...,
                   .drop = .drop)

  f.name <- as.character(substitute(.fun))

  comment(z) <- paste("Applied function: '", f.name, "'.\n", sep = "", comment(mspct))

  z
}

# generic_mspct methods -----------------------------------------------

#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#' @rdname  range.generic_spct
#'
range.generic_mspct <- function(..., na.rm = FALSE, idx = NULL) {
  mspct <- list(...)[[1]]
  if (is.null(idx)) {
    idx <- !is.null(names(mspct))
  }
  msdply(mspct = mspct, .fun = range, na.rm = na.rm, idx = idx)
}

#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#' @rdname  min.generic_spct
#'
min.generic_mspct <- function(..., na.rm = FALSE, idx = NULL) {
  mspct <- list(...)[[1]]
  if (is.null(idx)) {
    idx <- !is.null(names(mspct))
  }
  msdply(mspct = mspct, .fun = min, na.rm = na.rm, idx = idx)
}

#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#' @rdname  max.generic_spct
#'
max.generic_mspct <- function(..., na.rm = FALSE, idx = NULL) {
  mspct <- list(...)[[1]]
  if (is.null(idx)) {
    idx <- !is.null(names(mspct))
  }
  msdply(mspct = mspct, .fun = max, ..., na.rm = na.rm, idx = idx)
}

#' @describeIn stepsize  Method for "generic_mspct" objects.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
stepsize.generic_mspct <- function(x, ..., idx = !is.null(names(x))) {
  msdply(mspct = x, .fun = stepsize, ..., idx = idx)
}

#' @describeIn spread  Method for "generic_mspct" objects.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
spread.generic_mspct <- function(x, ..., idx = !is.null(names(x))) {
  msdply(mspct = x, .fun = spread, ..., idx = idx)
}

#' @describeIn midpoint Method for "generic_mspct" objects.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
midpoint.generic_mspct <- function(x, ..., idx = !is.null(names(x))) {
  msdply(mspct = x, .fun = midpoint, ..., idx = idx)
}

# when --------------------------------------------------------------------

#' @describeIn setWhenMeasured generic_mspct
#' @export
setWhenMeasured.generic_mspct <-
  function(x,
           when.measured = lubridate::now(),
           ...) {
    name <- substitute(x)
    stopifnot((lubridate::is.POSIXct(when.measured) && length(when.measured) == 1) ||
                is.list(when.measured))
    if (lubridate::is.POSIXct(when.measured) || length(when.measured) == 1) {
      if (is.list(when.measured)) {
        when.measured <- when.measured[[1]]
        stopifnot(lubridate::is.POSIXct(when.measured))
        lubridate::tz(when.measured) <- "UTC"
      }
      x <- msmsply(mspct = x, .fun = setWhenMeasured, when.measured = when.measured)
    } else if (length(when.measured) == length(x)) {
      for (i in 1:length(x)) {
        when <- when.measured[[i]]
        stopifnot(lubridate::is.POSIXct(when))
        lubridate::tz(when) <- "UTC"
        x[[i]] <- setWhenMeasured(x[[i]], when.measured = when)
      }
    }
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

#' @describeIn getWhenMeasured generic_mspct
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#' @note The method for collections of spectra returns the
#'   a data_frame with the correct times but they print in the local TZ as
#'   the "tzone" attribute is not supported by R's vectors of POSIXct objects
#'   or data frames.
#' @export
getWhenMeasured.generic_mspct <- function(x,
                                         ...,
                                         idx = !is.null(names(x))) {
  z <- msdply(mspct = x, .fun = getWhenMeasured, ..., idx = idx)
  attr(z, "tzone") <- "UTC"
  z
}

# where -------------------------------------------------------------------

#' @describeIn setWhereMeasured generic_mspct
#' @note Method for collections of spectra recycles the location information
#'   only if it is of length one.
#' @export
setWhereMeasured.generic_mspct <- function(x,
                                          where.measured = NA,
                                          lat = NA,
                                          lon = NA,
                                          ...) {
  name <- substitute(x)
  stopifnot(is.null(where.measured) || is.na(where.measured) ||
              is.data.frame(where.measured) ||
                (is.list(where.measured) && is.data.frame(where.measured[[1]])) )
  if (is.null(where.measured) ||
      (!is.na(where.measured) && is.data.frame(where.measured) && nrow(where.measured) == 1) ||
       (is.na(where.measured) && length(lat) == 1 && length(lon) == 1)) {
    x <- msmsply(mspct = x,
                 .fun = setWhereMeasured,
                 where.measured = where.measured,
                 lat = lat,
                 lon = lon)
  } else if (!is.na(where.measured) && !is.data.frame(where.measured) &&
             is.list(where.measured) && length(where.measured) == length(x)) {
    for (i in 1:length(x)) {
      x[[i]] <- setWhereMeasured(x[[i]], where.measured = where.measured[[i]])
    }
  } else if (!is.na(where.measured) && is.data.frame(where.measured) &&
                nrow(where.measured) == length(x)) {
    for (i in 1:length(x)) {
      x[[i]] <- setWhereMeasured(x[[i]], where.measured = where.measured[i, ])
    }
  } else if (is.na(where.measured) && length(lat) == length(x) && length(lon) == length(x)) {
    for (i in 1:length(x)) {
      x[[i]] <- setWhereMeasured(x[[i]], lon = lon[i], lat = lat[i])
    }
  } else {
    stop("Length of geocode information must be either 1, or equal to the number of spectra.")
  }
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' @describeIn getWhereMeasured generic_mspct
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#' @export
getWhereMeasured.generic_mspct <- function(x,
                                          ...,
                                          idx = !is.null(names(x))) {
  msdply(mspct = x, .fun = getWhereMeasured, ..., idx = idx)
}

# source_mspct methods -----------------------------------------------

#' @describeIn irrad  Calculates irradiance from a \code{source_mspct}
#'   object.
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
irrad.source_mspct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           allow.scaled = FALSE,
           ...,
           idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = irrad,
      w.band = w.band,
      unit.out = unit.out,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      allow.scaled = allow.scaled,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn q_irrad  Calculates photon (quantum) irradiance from a
#'   \code{source_mspct} object.
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
q_irrad.source_mspct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           allow.scaled = FALSE,
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = q_irrad,
      w.band = w.band,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      allow.scaled = allow.scaled,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn e_irrad  Calculates energy irradiance from a
#'   \code{source_mspct} object.
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
e_irrad.source_mspct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           allow.scaled = FALSE,
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = e_irrad,
      w.band = w.band,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      allow.scaled = allow.scaled,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn fluence Calculates fluence from a \code{source_mspct}
#'   object.
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
fluence.source_mspct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           exposure.time = NA,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           allow.scaled = FALSE,
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = fluence,
      w.band = w.band,
      unit.out = unit.out,
      exposure.time = exposure.time,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      allow.scaled = allow.scaled,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn e_fluence Calculates energy fluence from a \code{source_mspct}
#'   object.
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
e_fluence.source_mspct <-
  function(spct, w.band = NULL,
           exposure.time = NA,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           allow.scaled = FALSE,
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = e_fluence,
      w.band = w.band,
      exposure.time = exposure.time,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      allow.scaled = allow.scaled,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn q_fluence Calculates photon (quantum) fluence from a
#'   \code{source_mspct} object.
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
q_fluence.source_mspct <-
  function(spct, w.band = NULL,
           exposure.time = NA,
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           allow.scaled = FALSE,
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = q_fluence,
      w.band = w.band,
      exposure.time = exposure.time,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      allow.scaled = allow.scaled,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn q_ratio Calculates photon:photon from a \code{source_mspct}
#'   object.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
q_ratio.source_mspct <-
  function(spct,
           w.band.num = NULL, w.band.denom = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = q_ratio,
      w.band.num = w.band.num,
      w.band.denom = w.band.denom,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      idx = idx
    )
  }

#' @describeIn e_ratio Calculates energy:energy ratio from a \code{source_mspct}
#'   object.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
e_ratio.source_mspct <-
  function(spct,
           w.band.num = NULL, w.band.denom = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = e_ratio,
      w.band.num = w.band.num,
      w.band.denom = w.band.denom,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      idx = idx
    )
  }

#' @describeIn eq_ratio Calculates energy:photon from a \code{source_mspct}
#'   object.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
eq_ratio.source_mspct <-
  function(spct, w.band = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ...,
           idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = eq_ratio,
      w.band = w.band,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn qe_ratio Calculates photon:energy ratio from a
#'   \code{source_mspct} object.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
qe_ratio.source_mspct <-
  function(spct, w.band=NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges=getOption("photobiology.use.hinges", default = NULL),
           ...,
           idx = !is.null(names(spct))) {
    msdply(
      spct,
      .fun = qe_ratio,
      w.band = w.band,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#' @rdname color.source_spct
#'
color.source_mspct <- function(x, ..., idx = !is.null(names(x))) {
  msdply(mspct = x, color, ..., idx = idx)
}

# filter_mspct methods -----------------------------------------------

#' @describeIn transmittance Calculates transmittance from a \code{filter_mspct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
transmittance.filter_mspct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct)) ) {
    msdply(
      mspct = spct,
      .fun = transmittance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn absorptance Calculates absorptance from a \code{filter_mspct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
absorptance.filter_mspct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct)) ) {
    msdply(
      mspct = spct,
      .fun = absorptance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn absorbance Calculates absorbance from a \code{filter_mspct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
absorbance.filter_mspct <-
  function(spct, w.band=NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = absorbance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

# reflector_mspct methods -----------------------------------------------

#' @describeIn reflectance Calculates reflectance from a \code{reflector_mspct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
reflectance.reflector_mspct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = reflectance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

# object_mspct methods -----------------------------------------------

#' @describeIn transmittance Calculates transmittance from a \code{object_mspct}
#'
#' @export
#'
transmittance.object_mspct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct)) ) {
    msdply(
      mspct = spct,
      .fun = transmittance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }


#' @describeIn absorptance Calculates absorptance from a \code{object_mspct}
#'
#' @export
#'
absorptance.object_mspct <-
  function(spct, w.band=NULL,
           quantity="average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL),
           ..., idx = !is.null(names(spct)) ) {
    msdply(
      mspct = spct,
      .fun = absorptance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn reflectance Calculates reflectance from a \code{object_mspct}
#'
#' @export
#'
reflectance.object_mspct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges= getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = reflectance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn absorbance Calculates absorbance from a \code{object_mspct}
#'
#' @export
#'
absorbance.object_mspct <-
  function(spct, w.band=NULL,
           quantity="average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges=getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = absorbance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      col.names = names(w.band),
      idx = idx
    )
  }

# response_mspct methods -----------------------------------------------

#' @describeIn response Calculates response from a \code{response_mspct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
response.response_mspct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = response,
      w.band = w.band,
      unit.out = unit.out,
      quantity = quantity,
      time.unit = time.unit,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn q_response Calculates photon (quantum) response from a
#'   \code{response_mspct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
q_response.response_mspct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = q_response,
      w.band = w.band,
      quantity = quantity,
      time.unit = time.unit,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn e_response Calculates energy response from a
#'   \code{response_mspct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
e_response.response_mspct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = e_response,
      w.band = w.band,
      quantity = quantity,
      time.unit = time.unit,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }


#' Get the "mspct.version" attribute
#'
#' Funtion to read the "mspct.version" attribute of an existing generic_mspct
#' object.
#'
#' @param x a generic_mspct object
#'
#' @return numeric value
#'
#' @note if x is not a \code{generic_mspct} object, \code{NA} is returned,
#'   and if it the attribute is missing, zero is returned with a warning.
#'
#' @export
#'
getMspctVersion <- function(x) {
  if (is.any_mspct(x)) {
    version <- attr(x, "mspct.version", exact = TRUE)
    if (is.null(version)) {
      # need to handle objects created with old versions
      version <- 0L
    }
  } else {
    version <- NA
  }
  version
}

#' Check that the "mspct.version" attribute is set
#'
#' Funtion to check the "mspct.version" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_mspct object
#'
#' @return numeric value
#'
#' @note if x is not a \code{generic_mspct} object, \code{NA} is returned,
#'   and if it the attribute is missing, zero is returned with a warning.
#'
#' @keywords internal
#'
checkMspctVersion <- function(x) {
  version <- getMspctVersion(x)
  stopifnot(!is.na(version))
  if (version < 1L) {
    warning("The object '", as.character(substitute(x)),
            "' is corrupted")
  }
}


# print method ------------------------------------------------------------

#' Print a spectral collection object
#'
#' Print method for objects of the collection of spectra classes. Any object
#' of a class derived from \code{geenric_mspct} is printed with this  method.
#'
#' @param x generic_mspct A collection of spectra
#' @param ... not used in current version
#' @param n	Number of rows to show for each member spectrum.
#' @param width	Width of text output to generate.
#'
#' @return Returns \code{x} invisibly.
#'
#' @seealso \code{\link{print.generic_spct}} for details on \code{n}
#' and \code{width}.
#'
#' @export
#'
print.generic_mspct <- function(x, ..., n = NULL, width = NULL)  {
  cat("Object: ", class(x)[1], " ", dplyr::dim_desc(x), "\n", sep = "")
  member.names <- names(x)
  for (name in member.names) {
    cat("--- Member:", name, "---\n")
    print(x[[name]], n = n, width = width)
  }
  cat("--- END ---")
  invisible(x)
}

# convolute ---------------------------------------------------------------

#' Convolve function for collections of spectra
#'
#' Convolve function for collections of spectra which applies an operation on
#' all the individual members of the collection(s) of spectra.
#'
#' @param e1 an object of class \code{generic_mspct} or \code{generic_scpt} or
#'   \code{numeric}
#' @param e2 an object of class \code{generic_mspct} or \code{generic_scpt} or
#'   \code{numeric}
#' @param oper function, usually but not necesarily an operator with two
#'   arguments.
#' @param ... additional arguments passed to \code{oper} if present.
#'
#' @note At least one of e1 and e2 must be a \code{generic_mspct} object or
#'   derived.
#'
#' @export
#'
#' @family math operators and functions
#'
convolve_each <- function(e1, e2, oper = `*`, ...) {
  e3 <- list()
  if (is.any_mspct(e1) & !is.any_mspct(e2)) {
    for (spct.name in names(e1)) {
      e3[[spct.name]] <- oper(e1[[spct.name]], e2, ...)
    }
    z <- generic_mspct(e3, class = shared_member_class(e3),
                       ncol = ncol(e1),
                       byrow = attr(e1, "mspct.byrow", exact = TRUE))
  } else if (!is.any_mspct(e1) & is.any_mspct(e2)) {
    for (spct.name in names(e2)) {
      e3[[spct.name]] <- oper(e1, e2[[spct.name]], ...)
    }
    z <- generic_mspct(e3, class = shared_member_class(e3),
                       ncol = ncol(e2),
                       byrow = attr(e2, "mspct.byrow", exact = TRUE))
  } else if (is.any_mspct(e1) & is.any_mspct(e2)) {
    for (spct.name1 in names(e1)) {
      for (spct.name2 in names(e2)) {
        combined.name <- paste(spct.name1, spct.name2, sep = "_")
        e3[[combined.name]] <- oper(e1[[spct.name1]], e2[[spct.name2]], ...)
      }
     }
    z <- generic_mspct(e3, class = shared_member_class(e3),
                       ncol = nrow(e2),
                       byrow = FALSE)
    attr(z, "mspct.dimnames")  <- list(names(e1), names(e2))
  } else {
    stop("At least one of 'e1' and 'e2' should be a collection of spectra.")
  }
  z
}

# Extract ------------------------------------------------------------------

# $ operator for extraction does not need any wrapping as it always extracts
# single objects of the underlying classes (e.g. generic_spct)
# rather than collections of spectral objects.
#
# [ needs special handling as it can be used to extract members, or groups of
# members which must be returned as collections of spectral objects.
#
# In the case of replacement, collections of objects can easily become invalid,
# if the replacement or added member belongs to a class other than the expected
# one(s) for the collection.

#' Extract or replace members of a collection of spectra
#'
#' Just like extraction and replacement with indexes for base R lists, but
#' preserving the special attributes used in spectral classes.
#'
#' @param x	Collection of spectra object from which to extract member(s) or in
#'   which to replace member(s)
#' @param i Index specifying elements to extract or replace. Indices are numeric
#'   or character vectors. Please, see \code{\link[base]{Extract}} for
#'   more details.
#' @param drop If TRUE the result is coerced to the lowest possible dimension
#'   (see the examples). This only works for extracting elements, not for the
#'   replacement.
#'
#' @details This method is a wrapper on base R's extract method for lists that
#'   sets additional attributes used by these classes.
#'
#' @return An object of the same class as \code{x} but containing only the
#'   subset of members that are selected.
#'
#' @method [ generic_mspct
#' @export
#'
#' @rdname extract_mspct
#' @name Extract_mspct
#'
"[.generic_mspct" <-
  function(x, i, drop = NULL) {
    xx <- `[.listof`(x, i)
    generic_mspct(xx, class = class(x))
  }

# Not exported
# Check if class_spct is compatible with class_mspct
#
is.member_class <- function(l, x) {
  class(l)[1] == "generic_mscpt" && is.any_spct(x) ||
    sub("_mspct", "", class(l)[1], fixed = TRUE) == sub("_spct", "", class(x)[1], fixed = TRUE)
}

#' @param value	A suitable replacement value: it will be repeated a whole number
#'   of times if necessary and it may be coerced: see the Coercion section. If
#'   NULL, deletes the column if a single column is selected.
#'
#' @export
#' @method [<- generic_mspct
#' @rdname extract_mspct
#'
"[<-.generic_mspct" <- function(x, i, value) {
  # could be improved to accept derived classes as valid for replacement.
  stopifnot(class(x) == class(value))
  # could not find a better way of avoiding infinite recursion as '[[<-' is
  # a primitive with no explicit default method.
  old.class <- class(x)
  class(x) <- "list"
  x[i] <- value
  class(x) <- old.class
  x
}

#' @param name A literal character string or a name (possibly backtick quoted).
#'   For extraction, this is normally (see under ‘Environments’) partially
#'   matched to the names of the object.
#'
#' @export
#' @method $<- generic_mspct
#' @rdname extract_mspct
#'
"$<-.generic_mspct" <- function(x, name, value) {
  x[[name]] <- value
}

#' @export
#' @method [[<- generic_mspct
#' @rdname extract_mspct
#'
"[[<-.generic_mspct" <- function(x, name, value) {
  stopifnot(is.member_class(x, value) || is.null(value))
  # could not find a better way of avoiding infinite recursion as '[[<-' is
  # a primitive with no explicit default method.
  if (is.character(name) && !(name %in% names(x)) ) {
    if (ncol(x) == 1) {
      dimension <- c(nrow(x) + 1, 1)
  } else {
      stop("Appending to a matrix-like collection not supported.")
    }
  } else if (is.numeric(name) && (name > length(x)) ) {
    stop("Appending to a collection using numeric indexing not supported.")
  } else if (is.null(value)) {
    if (ncol(x) != 1) {
      stop("Deleting members from a matrix-like collection not supported.")
    } else {
      dimension <- attr(x, "mspct.dim", exact = TRUE)
    }
  } else {
    dimension <- attr(x, "mspct.dim", exact = TRUE)
  }
  old.class <- class(x)
  old.byrow <- attr(x, "mspct.byrow", exact = TRUE)
  class(x) <- "list"
  x[[name]] <- value
  class(x) <- old.class
  attr(x, "mspct.dim") <- dimension
  attr(x, "mspct.byrow") <- old.byrow
  x
}

# Combine -----------------------------------------------------------------

#' Combine collections of spectra
#'
#' Combine two or more generic_mspct objects into a single object.
#'
#' @param ... one or more generic_mspct objects to combine.
#' @param recursive logical ignored as nesting of collections of spectra is
#' not supported.
#' @param ncol numeric Virtual number of columns
#' @param byrow logical When object has two dimensions, how to map member
#' objects to columns and rows.
#'
#' @return A collection of spectra object belonging to the most derived class
#' shared among the combined objects.
#'
#' @export
#' @method c generic_mspct
#'
c.generic_mspct <- function(..., recursive = FALSE, ncol = 1, byrow = FALSE) {
  l <- list(...)
  shared.class <- shared_member_class(l, target.set = mspct_classes())
  stopifnot(length(shared.class) > 0)
  shared.class <- shared.class[1]
  ul <- unlist(l, recursive = FALSE)
  do.call(shared.class, list(l = ul, ncol = ncol, byrow = byrow))
}

