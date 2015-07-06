#' Multi-sppct summary methods
#'
#' Functions
#'
#' @param mspct an object of class generic_multi_spct or a derived class
#' @param f a function
#' @param ... other arguments passed to irrad.source.spct
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
f_multi_spct <- function(mspct, f, ..., idx = !is.null(names(mspct))) {
  z <- NULL

  z <- lapply(mspct, f, ...)
  namesz <- names(z[[1]])
  z <- unlist(z, recursive = FALSE, use.names = FALSE)

  nspct <- length(mspct)
  nz <- length(z)
  nqty <- nz %/% nspct
  ncol <- attr(mspct, "ncol", exact = TRUE)
  nrow <- nspct / ncol

  stopifnot(nspct %% ncol == 0)

  rows <- rep(1:nrow, ncol)
  cols <- rep(1:ncol, rep(nrow, ncol))

  if (nz / nspct > 1) {
    z <- matrix(z, ncol = nqty, byrow = TRUE)
  }
  df <- data.frame(row = rows, col = cols, z = z)
  if (idx) {
    df[["idx"]] <- names(mspct)
  }

  qty.names <- paste(as.character(substitute(f)), gsub(" ", "", namesz), sep = "_")

  setnames(df, (1:nqty) + 2, qty.names)
  df
}

# generic_multi_spct methods -----------------------------------------------

#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#' @rdname  range.generic_spct
#'
range.generic_multi_spct <- function(..., na.rm = FALSE, idx = !is.null(names(spct))) {
  mspct <- c(...)
  f_multi_spct(mspct = mspct, f = range, na.rm = na.rm, idx = idx)
}

#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#' @rdname  min.generic_spct
#'
min.generic_multi_spct <- function(..., na.rm = FALSE, idx = !is.null(names(spct))) {
  mspct <- c(...)
  f_multi_spct(mspct = mspct, min, na.rm = na.rm, idx = idx)
}

#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#' @rdname  max.generic_spct
#'
max.generic_multi_spct <- function(..., na.rm = FALSE, idx = !is.null(names(spct))) {
  mspct <- c(...)
  f_multi_spct(mspct = mspct, max, ..., na.rm = na.rm, idx = idx)
}

#' @describeIn stepsize  Method for "generic_multi_spct" objects for generic function.
#'
##' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
stepsize.generic_multi_spct <- function(x, ..., idx = !is.null(names(spct))) {
  f_multi_spct(mspct = x, stepsize, ..., idx = idx)
}

#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#' @rdname  labels.generic_spct
#'
labels.generic_multi_spct <- function(object, ..., idx = !is.null(names(spct))) {
  f_multi_spct(object, labels, ..., idx = idx)
}

# source_multi_spct methods -----------------------------------------------

#' @describeIn irrad  Calculates irradiance from a \code{source_multi_spct}
#'   object.
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
irrad.source_multi_spct <-
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
  f_multi_spct(mspct = spct, f = irrad,
               w.band = w.band, unit.out = unit.out,
               wb.trim = wb.trim, use.cached.mult = use.cached.mult,
               use.hinges = use.hinges, allow.scaled = allow.scaled, idx = idx)
}

#' @describeIn q_irrad  Calculates photon (quantum) irradiance from a
#'   \code{source_multi_spct} object.
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
q_irrad.source_multi_spct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           allow.scaled = FALSE,
           ..., idx = !is.null(names(spct))) {
  f_multi_spct(mspct = spct, f = q_irrad,
               w.band = w.band,
               wb.trim = wb.trim, use.cached.mult = use.cached.mult,
               use.hinges = use.hinges, allow.scaled = allow.scaled, idx = idx)
}

#' @describeIn e_irrad  Calculates energy irradiance from a
#'   \code{source_multi_spct} object.
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
e_irrad.source_multi_spct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           allow.scaled = FALSE,
           ..., idx = !is.null(names(spct))) {
  f_multi_spct(mspct = spct, f = e_irrad,
               w.band = w.band,
               wb.trim = wb.trim, use.cached.mult = use.cached.mult,
               use.hinges = use.hinges, allow.scaled = allow.scaled, idx = idx)
}

#' @describeIn fluence Calculates fluence from a \code{source_multi_spct}
#'   object.
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
fluence.source_multi_spct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           exposure.time = NA,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           allow.scaled = FALSE,
           ..., idx = !is.null(names(spct))) {
    f_multi_spct(mspct = spct, f = fluence,
                 w.band = w.band, unit.out = unit.out,
                 exposure.time = exposure.time,
                 wb.trim = wb.trim, use.cached.mult = use.cached.mult,
                 use.hinges = use.hinges, allow.scaled = allow.scaled, idx = idx)
}

#' @describeIn e_fluence Calculates energy fluence from a \code{source_multi_spct}
#'   object.
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
e_fluence.source_multi_spct <-
  function(spct, w.band = NULL,
           exposure.time = NA,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           allow.scaled = FALSE,
           ..., idx = !is.null(names(spct))) {
    f_multi_spct(mspct = spct, f = e_fluence,
                 w.band = w.band,
                 exposure.time = exposure.time,
                 wb.trim = wb.trim, use.cached.mult = use.cached.mult,
                 use.hinges = use.hinges, allow.scaled = allow.scaled, idx = idx)
  }

#' @describeIn q_fluence Calculates photon (quantum) fluence from a
#'   \code{source_multi_spct} object.
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
q_fluence.source_multi_spct <-
  function(spct, w.band = NULL,
           exposure.time = NA,
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           allow.scaled = FALSE,
           ..., idx = !is.null(names(spct))) {
    f_multi_spct(mspct = spct, f = q_fluence,
                 w.band = w.band,
                 exposure.time = exposure.time,
                 wb.trim = wb.trim, use.cached.mult = use.cached.mult,
                 use.hinges = use.hinges, allow.scaled = allow.scaled, idx = idx)
  }

#' @describeIn q_ratio Calculates photon:photon from a \code{source_multi_spct}
#'   object.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
q_ratio.source_multi_spct <-
  function(spct,
           w.band.num = NULL, w.band.denom = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
  f_multi_spct(mspct = spct, f = q_ratio,
               w.band.num = w.band.num, w.band.denom = w.band.denom,
               wb.trim = wb.trim, use.cached.mult = use.cached.mult,
               use.hinges = use.hinges, idx = idx)
}

#' @describeIn e_ratio Calculates energy:energy ratio from a \code{source_multi_spct}
#'   object.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
e_ratio.source_multi_spct <-
  function(spct,
           w.band.num = NULL, w.band.denom = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
  f_multi_spct(mspct = spct, f = e_ratio,
               w.band.num = w.band.num, w.band.denom = w.band.denom,
               wb.trim = wb.trim, use.cached.mult = use.cached.mult,
               use.hinges = use.hinges, idx = idx)
}

#' @describeIn eq_ratio Calculates energy:photon from a \code{source_multi_spct}
#'   object.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
eq_ratio.source_multi_spct <-
  function(spct, w.band = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ...,
           idx = !is.null(names(spct))) {
  f_multi_spct(mspct = spct, f = eq_ratio,
               w.band = w.band,
               wb.trim = wb.trim, use.cached.mult = use.cached.mult,
               use.hinges = use.hinges, idx = idx)
}

#' @describeIn qe_ratio Calculates photon:energy ratio from a
#'   \code{source_multi_spct} object.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
qe_ratio.source_multi_spct <-
  function(spct, w.band=NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges=getOption("photobiology.use.hinges", default = NULL),
           ...,
           idx = !is.null(names(spct))) {
  f_multi_spct(spct, qe_ratio,
               w.band = w.band,
               wb.trim = wb.trim, use.cached.mult = use.cached.mult,
               use.hinges = use.hinges, idx = idx)
}

#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#' @rdname color.source_spct
#'
color.source_multi_spct <- function(x, ..., idx = !is.null(names(spct))) {
  f_multi_spct(mspct = x, color, ..., idx = idx)
}

# filter_multi_spct methods -----------------------------------------------

#' @describeIn transmittance Calculates transmittance from a
#'   \code{filter_multi_spct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
transmittance.filter_multi_spct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct)) ) {
    f_multi_spct(mspct = spct, f = transmittance,
                 w.band = w.band,
                 quantity = quantity,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
  }

#' @describeIn absorptance Calculates absorptance from a
#'   \code{filter_multi_spct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
absorptance.filter_multi_spct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct)) ) {
    f_multi_spct(mspct = spct, f = absorptance,
                 w.band = w.band,
                 quantity = quantity,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
  }

#' @describeIn absorbance Calculates absorbance from a \code{filter_multi_spct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
absorbance.filter_multi_spct <-
  function(spct, w.band=NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
  f_multi_spct(mspct = spct, f = absorbance,
               w.band = w.band,
               quantity = quantity,
               wb.trim = wb.trim,
               use.hinges = use.hinges,
               idx = idx)
}

# reflector_multi_spct methods -----------------------------------------------

#' @describeIn reflectance Calculates reflectance from a
#'   \code{reflector_multi_spct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
reflectance.reflector_multi_spct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
  f_multi_spct(mspct = spct, f = reflectance,
               w.band = w.band,
               quantity = quantity,
               wb.trim = wb.trim,
               use.hinges = use.hinges,
               idx = idx)
}

# object_multi_spct methods -----------------------------------------------

#' @describeIn transmittance Calculates transmittance from a
#'   \code{object_multi_spct}
#'
#' @export
#'
transmittance.object_multi_spct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct)) ) {
    f_multi_spct(mspct = spct, f = transmittance,
                 w.band = w.band,
                 quantity = quantity,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
  }


#' @describeIn absorptance Calculates absorptance from a
#'   \code{object_multi_spct}
#'
#' @export
#'
absorptance.object_multi_spct <-
  function(spct, w.band=NULL,
           quantity="average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL),
           ..., idx = !is.null(names(spct)) ) {
    f_multi_spct(mspct = spct, f = absorptance,
                 w.band = w.band,
                 quantity = quantity,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
  }

#' @describeIn reflectance Calculates reflectance from a
#'   \code{object_multi_spct}
#'
#' @export
#'
reflectance.object_multi_spct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges= getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    f_multi_spct(mspct = spct, f = reflectance,
                 w.band = w.band,
                 quantity = quantity,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
  }

#' @describeIn absorbance Calculates absorbance from a \code{object_multi_spct}
#'
#' @export
#'
absorbance.object_multi_spct <-
  function(spct, w.band=NULL,
           quantity="average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL),
           ..., idx = !is.null(names(spct))) {
    f_multi_spct(mspct = spct, f = absorbance,
                 w.band = w.band,
                 quantity = quantity,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
  }

# response_multi_spct methods -----------------------------------------------

#' @describeIn response Calculates response from a \code{response_multi_spct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
response.response_multi_spct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
  f_multi_spct(mspct = spct, f = response,
               w.band = w.band, unit.out = unit.out,
               quantity = quantity,
               time.unit = time.unit,
               wb.trim = wb.trim,
               use.hinges = use.hinges,
               idx = idx)
}

#' @describeIn q_response Calculates photon (quantum) response from a
#'   \code{response_multi_spct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
q_response.response_multi_spct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
  f_multi_spct(mspct = spct, f = q_response,
               w.band = w.band,
               quantity = quantity,
               time.unit = time.unit,
               wb.trim = wb.trim,
               use.hinges = use.hinges,
               idx = idx)
}

#' @describeIn e_response Calculates energy response from a
#'   \code{response_multi_spct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
e_response.response_multi_spct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    f_multi_spct(mspct = spct, f = e_response,
                 w.band = w.band,
                 quantity = quantity,
                 time.unit = time.unit,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
  }

