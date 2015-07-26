#' Multi-spct transform methods
#'
#' Apply a function returning an object of the same class as its first argument.
#'
#' @param mspct an object of class generic_mspct or a derived class
#' @param f a function
#' @param ... other arguments passed to f
#'
#' @return an object of the same class as mspct
#'
#' @export
#'
mutate_mspct <- function(mspct, f, ...) {
  stopifnot(is.any_mspct(mspct))
  mspct.class <- class(mspct)[1]

  y <- llply(mspct, f, ...)

  stopifnot(length(y) == length(mspct))

  generic_mspct(l = y,
                class = mspct.class,
                byrow = attr(mspct, "byrow", exact = TRUE),
                ncol = ncol(mspct))
}

#' Multi-spct summary methods
#'
#' Functions
#'
#' @param mspct an object of class generic_mspct or a derived class
#' @param f a function
#' @param ... other arguments passed to f
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @return a data frame
#'
#' @export
#'
f_mspct <- function(mspct, f, ..., idx = NULL) {
  stopifnot(is.any_mspct(mspct))

  if ( (is.logical(idx) && idx) ||
       (is.null(idx) && !any(is.null(names(mspct)))) ) {
    idx <- "spct.idx"
  }

  z <- plyr::ldply(.data = mspct,
                   .fun = f,
                   ...,
                   .id = idx )

  f.name <- as.character(substitute(f))

  if (f.name %in% c("min", "max", "range", "spread", "midpoint", "stepsize")) {
    qty.names <- switch(f.name,
                        min = "min.wl",
                        max = "max.wl",
                        range = c("min.wl", "max.wl"),
                        spread = "spread.wl",
                        midpoint = "midpoint.wl",
                        stepsize = c("min.step.wl", "max.step.wl") )
  } else {
    qty.names <- paste(f.name,
                       gsub(" ", "", names(z)[-1]),
                       sep = "_")
  }

  names(z) <- c(names(z)[1], qty.names)

  mspct.nrow <- nrow(mspct)
  mspct.ncol <- ncol(mspct)

  if (mspct.ncol > 1) {
    if (byrow) {
      z$col <- rep(1:mspct.ncol, mspct.nrow)
      z$row <- rep(1:mspct.nrow, rep(mspct.col, mspct.nrow))
    } else {
      z$col <- rep(1:mspct.ncol, rep(mspct.nrow, mspct.ncol))
      z$row <- rep(1:mspct.nrow, mspct.col)
    }
  }
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
  f_mspct(mspct = mspct, f = range, na.rm = na.rm, idx = idx)
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
  f_mspct(mspct = mspct, f = min, na.rm = na.rm, idx = idx)
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
  f_mspct(mspct = mspct, f = max, ..., na.rm = na.rm, idx = idx)
}

#' @describeIn stepsize  Method for "generic_mspct" objects.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
stepsize.generic_mspct <- function(x, ..., idx = !is.null(names(x))) {
  f_mspct(mspct = x, f = stepsize, ..., idx = idx)
}

#' @describeIn spread  Method for "generic_mspct" objects.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
spread.generic_mspct <- function(x, ..., idx = !is.null(names(x))) {
  f_mspct(mspct = x, f = spread, ..., idx = idx)
}

#' @describeIn midpoint Method for "generic_mspct" objects.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
midpoint.generic_mspct <- function(x, ..., idx = !is.null(names(x))) {
  f_mspct(mspct = x, f = midpoint, ..., idx = idx)
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
    f_mspct(mspct = spct, f = irrad,
            w.band = w.band, unit.out = unit.out,
            wb.trim = wb.trim, use.cached.mult = use.cached.mult,
            use.hinges = use.hinges, allow.scaled = allow.scaled, idx = idx)
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
  f_mspct(mspct = spct, f = q_irrad,
               w.band = w.band,
               wb.trim = wb.trim, use.cached.mult = use.cached.mult,
               use.hinges = use.hinges, allow.scaled = allow.scaled, idx = idx)
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
  f_mspct(mspct = spct, f = e_irrad,
               w.band = w.band,
               wb.trim = wb.trim, use.cached.mult = use.cached.mult,
               use.hinges = use.hinges, allow.scaled = allow.scaled, idx = idx)
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
    f_mspct(mspct = spct, f = fluence,
                 w.band = w.band, unit.out = unit.out,
                 exposure.time = exposure.time,
                 wb.trim = wb.trim, use.cached.mult = use.cached.mult,
                 use.hinges = use.hinges, allow.scaled = allow.scaled, idx = idx)
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
    f_mspct(mspct = spct, f = e_fluence,
                 w.band = w.band,
                 exposure.time = exposure.time,
                 wb.trim = wb.trim, use.cached.mult = use.cached.mult,
                 use.hinges = use.hinges, allow.scaled = allow.scaled, idx = idx)
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
    f_mspct(mspct = spct, f = q_fluence,
                 w.band = w.band,
                 exposure.time = exposure.time,
                 wb.trim = wb.trim, use.cached.mult = use.cached.mult,
                 use.hinges = use.hinges, allow.scaled = allow.scaled, idx = idx)
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
  f_mspct(mspct = spct, f = q_ratio,
               w.band.num = w.band.num, w.band.denom = w.band.denom,
               wb.trim = wb.trim, use.cached.mult = use.cached.mult,
               use.hinges = use.hinges, idx = idx)
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
  f_mspct(mspct = spct, f = e_ratio,
               w.band.num = w.band.num, w.band.denom = w.band.denom,
               wb.trim = wb.trim, use.cached.mult = use.cached.mult,
               use.hinges = use.hinges, idx = idx)
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
  f_mspct(mspct = spct, f = eq_ratio,
               w.band = w.band,
               wb.trim = wb.trim, use.cached.mult = use.cached.mult,
               use.hinges = use.hinges, idx = idx)
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
  f_mspct(spct, qe_ratio,
               w.band = w.band,
               wb.trim = wb.trim, use.cached.mult = use.cached.mult,
               use.hinges = use.hinges, idx = idx)
}

#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#' @rdname color.source_spct
#'
color.source_mspct <- function(x, ..., idx = !is.null(names(spct))) {
  f_mspct(mspct = x, color, ..., idx = idx)
}

# filter_mspct methods -----------------------------------------------

#' @describeIn transmittance Calculates transmittance from a
#'   \code{filter_mspct}
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
    f_mspct(mspct = spct, f = transmittance,
                 w.band = w.band,
                 quantity = quantity,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
  }

#' @describeIn absorptance Calculates absorptance from a
#'   \code{filter_mspct}
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
    f_mspct(mspct = spct, f = absorptance,
                 w.band = w.band,
                 quantity = quantity,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
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
  f_mspct(mspct = spct, f = absorbance,
               w.band = w.band,
               quantity = quantity,
               wb.trim = wb.trim,
               use.hinges = use.hinges,
               idx = idx)
}

# reflector_mspct methods -----------------------------------------------

#' @describeIn reflectance Calculates reflectance from a
#'   \code{reflector_mspct}
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
  f_mspct(mspct = spct, f = reflectance,
               w.band = w.band,
               quantity = quantity,
               wb.trim = wb.trim,
               use.hinges = use.hinges,
               idx = idx)
}

# object_mspct methods -----------------------------------------------

#' @describeIn transmittance Calculates transmittance from a
#'   \code{object_mspct}
#'
#' @export
#'
transmittance.object_mspct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct)) ) {
    f_mspct(mspct = spct, f = transmittance,
                 w.band = w.band,
                 quantity = quantity,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
  }


#' @describeIn absorptance Calculates absorptance from a
#'   \code{object_mspct}
#'
#' @export
#'
absorptance.object_mspct <-
  function(spct, w.band=NULL,
           quantity="average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL),
           ..., idx = !is.null(names(spct)) ) {
    f_mspct(mspct = spct, f = absorptance,
                 w.band = w.band,
                 quantity = quantity,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
  }

#' @describeIn reflectance Calculates reflectance from a
#'   \code{object_mspct}
#'
#' @export
#'
reflectance.object_mspct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges= getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    f_mspct(mspct = spct, f = reflectance,
                 w.band = w.band,
                 quantity = quantity,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
  }

#' @describeIn absorbance Calculates absorbance from a \code{object_mspct}
#'
#' @export
#'
absorbance.object_mspct <-
  function(spct, w.band=NULL,
           quantity="average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL),
           ..., idx = !is.null(names(spct))) {
    f_mspct(mspct = spct, f = absorbance,
                 w.band = w.band,
                 quantity = quantity,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
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
  f_mspct(mspct = spct, f = response,
               w.band = w.band, unit.out = unit.out,
               quantity = quantity,
               time.unit = time.unit,
               wb.trim = wb.trim,
               use.hinges = use.hinges,
               idx = idx)
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
  f_mspct(mspct = spct, f = q_response,
               w.band = w.band,
               quantity = quantity,
               time.unit = time.unit,
               wb.trim = wb.trim,
               use.hinges = use.hinges,
               idx = idx)
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
    f_mspct(mspct = spct, f = e_response,
                 w.band = w.band,
                 quantity = quantity,
                 time.unit = time.unit,
                 wb.trim = wb.trim,
                 use.hinges = use.hinges,
                 idx = idx)
  }

#' Get the "mspct.version" attribute
#'
#' Funtion to read the "mspct.version" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_mspct object
#'
#' @return numeric value
#'
#' @note if x is not a \code{generic_spct} object, \code{NA} is returned,
#'   and if it the attribute is missing, zero is returned with a warning.
#'
#' @export
#'
getMspctVersion <- function(x) {
  if (is.any_spct(x)) {
    version <- attr(x, "spct.version", exact = TRUE)
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

