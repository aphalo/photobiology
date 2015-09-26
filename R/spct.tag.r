
# tag ---------------------------------------------------------------------

#' Tag a spectrum
#'
#' Spectra are tagged by adding variables and attributes containing color
#' definitions, labels, and a factor following the wavebands given in
#' \code{w.band}.
#'
#' @param x an R object
#' @param ... not used in current version
#'
#' @export tag
#'
#' @family tagging and related functions
#'
tag <- function(x, ...) UseMethod("tag")

#' @describeIn tag Default method for generic
#'
#' @export
#'
tag.default <- function(x, ...) {
  warning("'tag' is not defined for objects of class ", class(x)[1])
  return(x)
}

#' @describeIn tag Tag one of \code{generic_spct}, and derived classes including
#'   \code{source_spct}, \code{filter_spct}, \code{reflector_spct},
#'   \code{object_spct}, and \code{rspponse.spct}.
#'
#' @param w.band waveband or list of waveband objects The waveband(s) determine
#'   the region(s) of the spectrum that are tagged
#' @param wb.trim logical Flag telling if wavebands crossing spectral data
#'   boundaries are trimmed or ignored
#' @param use.hinges logical Flag indicating whether to use hinges to reduce
#'   interpolation errors
#' @param short.names logical Flag indicating whether to use short or long names
#'   for wavebands
#' @param byref logical Flag indicating if new object will be created \emph{by
#'   reference} or \emph{by copy} of \code{x}
#' @export
#'
#' @note \code{NULL} as \code{w.band} argument does not add any new tags,
#'   instead it removes existing tags if present. \code{NA}, the default, as
#'   \code{w.band} argument removes existing waveband tags if present and
#'   sets the \code{wl.color} variable. If a waveband object or a list of
#'   wavebands is supplied as argument then tagging is based on them, and
#'   \code{wl.color} is also set.
#'
#' @examples
#'
#' tag(sun.spct)
#' tag(sun.spct, list(A = waveband(c(300,3005))))
#'
tag.generic_spct <- function(x,
                             w.band = NA,
                             wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
                             use.hinges = TRUE,
                             short.names = TRUE,
                             byref = FALSE, ...) {
  name <- substitute(x)
  if (is_tagged(x)) {
    warning("Overwriting old tags in spectrum")
    untag(x, byref = TRUE)
  }
 # we add a waveband for the whole spectrum
  if (is.null(w.band)) {
     return(x)
  }
  if (!is.null(w.band) && is.na(w.band[1])) {
    x[["wl.color"]] <- w_length2rgb(x[["w.length"]])
    tag.data <- list(wl.color = TRUE)
    attr(x, "spct.tags") <- tag.data
    return(x)
  }
  if (!is.null(w.band) && is(w.band, "waveband")) {
    # if the argument is a single w.band, we enclose it in a list
    # so that the for loop works as expected.This is a bit of a
    # cludge but lets us avoid treating it as a special case
    w.band <- list(w.band)
  }
  # we delete or trim the wavebands that are not fully within the
  # spectral data wavelngth range
  w.band <- trim_waveband(w.band = w.band, range = x, trim = wb.trim)
  # we check if the list elements are named, if not we set a flag
  # and an empty vector that will be later filled in with data from
  # the waveband definitions.
  wbs.number <- length(w.band) # number of wavebands in list
  wbs.name <- names(w.band) # their names in the list
  if (is.null(wbs.name)) {
    wbs.name <- character(wbs.number)
  }
  # if the w.band includes 'hinges' we insert them
  # choose whether to use hinges or not
  # if the user has specified its value, we leave it alone
  # but if it was not requested, we decide whether to use
  # it or not based of the wavelength resolution of the
  # spectrum. This will produce small errors for high
  # spectral resolution data, and speed up the calculations
  # a lot in such cases
  if (is.null(use.hinges)) {
    use.hinges <- auto_hinges(x)
  }
  # we collect all hinges and insert them in one go
  # this may alter a little the returned values
  # but should be faster
  if (use.hinges) {
    all.hinges <- NULL
    for (wb in w.band) {
      if (!is.null(wb$hinges) & length(wb$hinges) > 0) {
        all.hinges <- c(all.hinges, wb[["hinges"]])
      }
    }
    if (!is.null(all.hinges)) {
      x <- insert_spct_hinges(x, all.hinges)
    }
  }

  # We iterate through the list of wavebands adding the tags
  wbs.rgb <- character(wbs.number)
  wbs.wl.low <- wbs.wl.high <- numeric(wbs.number)
  i <- 0L
  for (wb in w.band) {
    i <- i + 1L
    if (wbs.name[i] == "") {
      if (short.names) {
        name.temp <- labels(wb)[["label"]]
        wbs.name[i] <- ifelse(grepl("^range.", name.temp, ignore.case = TRUE),
                              paste("wb", i, sep = ""),
                              name.temp)
      } else {
        wbs.name[i] <- labels(wb)[["name"]]
      }
    }
    wbs.wl.low[i] <- min(wb)
    wbs.wl.high[i] <- max(wb)
    wbs.rgb[i] <- color(wb)[1]
  }
  n <- i
  x[["wl.color"]] <- NA
  x[["wb.f"]] <- NA
  for (i in 1L:n) {
    selector <- x[["w.length"]] >= wbs.wl.low[i] & x[["w.length"]] < wbs.wl.high[i]
    x[selector, "wb.f"] <- wbs.name[i]
  }
  x[["wl.color"]] <-  w_length2rgb(x[["w.length"]])
  x[["wb.f"]] <- factor(x[["wb.f"]], levels = wbs.name)
  tag.data <- list(time.unit = getTimeUnit(x),
                   wb.key.name = "Bands",
                   wl.color = TRUE,
                   wb.color = TRUE,
                   wb.num = n,
                   wb.colors = wbs.rgb[1:n],
                   wb.names = wbs.name[1:n],
                   wb.list = w.band)
  attr(x, "spct.tags") <- tag.data
  # to work by reference we need to assign the new data frame to the old one
  if (byref & is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}


# wavebands -> tagged spectrum ----------------------------------------------

#' Create spectrum from wavebands
#'
#' Create a generic_spct object with wavelengths from wavebands in a list.
#'
#' @param w.band waveband or list of waveband objects The waveband(s) determine
#'   the wavelengths in variable \code{w.length} of the returned spectrum
#' @export
#'
#' @return A generic.spectrum object, with columns w.length set to the
#'   \emph{union} of all  boundaries and hinges defined in the waveband(s).
#'   Different spectral data variables are set to zero and added making the
#'   returned value compatible with classes derived from \code{generic_spct}.
#'
#' @family tagging and related functions
#'
wb2spct <- function(w.band) {
  if (is.waveband(w.band)) {
    w.band <- list(w.band)
  }
  w.length <- numeric(0)
  for (wb in w.band) {
    stopifnot(is.waveband(wb))
    w.length <- c(w.length, wb$hinges)
  }
  if (is.null(w.length) || length(w.length) < 2) {
    new.spct <- dplyr::data_frame(w.length = numeric(0), s.e.irrad = numeric(0), s.q.irrad = numeric(0),
                                  Tfr = numeric(0), Rfl = numeric(0), s.e.response = numeric(0))
  } else {
    w.length <- unique(sort(w.length))
    new.spct <- dplyr::data_frame(w.length = w.length, s.e.irrad = 0, s.q.irrad = 0,
                                  Tfr = 0, Rfl = 0, s.e.response = 0)
  }
  setGenericSpct(new.spct)
  new.spct
}

#' Create tagged spectrum from wavebands
#'
#' Create a tagged \code{generic_spct} object with wavelengths from the range of
#' wavebands in a list, and names of the same bands as factor levels, and
#' corresponding color definitions. The spectrum is not suitable for plotting
#' labels, symbols, rectangles or similar, as the midpoint of each waveband is
#' not added to the spectrum.
#'
#' @param w.band waveband or list of waveband objects The waveband(s) determine
#'   the region(s) of the spectrum that are tagged and the wavelengths returned
#'   in variable \code{w.length}
#' @param use.hinges logical Flag indicating whether to use hinges to reduce
#'   interpolation errors
#' @param short.names logical Flag indicating whether to use short or long names
#'   for wavebands
#' @param ... not used in current version
#' @export
#'
#' @return A spectrum as returned by \code{\link{wb2spct}} but additionally
#'   tagged using function \code{\link{tag}}
#'
#' @family tagging and related functions
#'
wb2tagged_spct <-
  function(w.band, use.hinges = TRUE, short.names = TRUE, ...) {
  new.spct <- wb2spct(w.band)
  tag(new.spct, w.band, use.hinges, short.names, byref = TRUE)
  new.spct[["y"]] <- 0
  return(new.spct)
}

#' Create tagged spectrum from wavebands
#'
#' Create a generic_spct object with wavelengths from the range of wavebands in
#' a list. The spectrum is suitable for plotting labels, symbols, rectangles or
#' similar, as the midpoint of each waveband is added to the spectrum.
#'
#' @param w.band waveband or list of waveband objects The waveband(s) determine
#'   the wavelengths in variable \code{w.length} of the returned spectrum
#' @param short.names logical Flag indicating whether to use short or long names
#'   for wavebands
#' @export
#'
#' @return A \code{generic.spectrum} object, with columns w.length, wl.low,
#'   wl.hi, and wl.color. The w.length values are such that at the midppoint of
#'   the wavebands, wl.low and wl.high give the boundaries of the wavebands, and
#'   wl.color the color definition corresponding to the wavebands. Different
#'   spectral data variables are set to zero and added making the returned value
#'   compatible with classes derived from \code{generic_spct}.
#'
#' @family tagging and related functions
#'
wb2rect_spct <- function(w.band, short.names = TRUE) {
  if (is(w.band, "waveband")) {
    w.band <- list(w.band)
  }
  wbs.number <- length(w.band) # number of wavebands in list
  wbs.name <- names(w.band)
  if (is.null(wbs.name)) {
    wbs.name <- character(wbs.number)
  }
  wbs.wl.mid <- wbs.wl.high <- wbs.wl.low <- numeric(wbs.number)
  wbs.rgb <- character(wbs.number)
  i <- 0L
  for (wb in w.band) {
    i <- i + 1L
    if (wbs.name[i] == "") {
      if (short.names) {
        name.temp <- labels(wb)[["label"]]
        wbs.name[i] <- ifelse(grepl("^range.", name.temp, ignore.case = TRUE),
                              paste("wb", i, sep = ""),
                              name.temp)
      } else {
        wbs.name[i] <- labels(wb)[["name"]]
      }
    }
    wbs.wl.low[i] <- min(wb)
    wbs.wl.mid[i] <- midpoint(wb)
    wbs.wl.high[i] <- max(wb)
    wbs.rgb[i] <- color(wb)[1]
  }
  new.spct <- dplyr::data_frame(w.length = wbs.wl.mid,
                                s.e.irrad = 0, s.q.irrad = 0, Tfr = 0, Rfl = 0, s.e.response = 0,
                                wl.color = w_length2rgb(wbs.wl.mid),
                                wb.f = factor(wbs.name, levels = wbs.name),
                                wl.high = wbs.wl.high, wl.low = wbs.wl.low,
                                y = 0)
  setGenericSpct(new.spct)
  tag.data <- list(time.unit = "none",
                   wb.key.name = "Bands",
                   wl.color = TRUE,
                   wb.color = TRUE,
                   wb.num = wbs.number,
                   wb.colors = wbs.rgb,
                   wb.names = wbs.name,
                   wb.list = w.band)
  attr(new.spct, "spct.tags") <- tag.data

  return(new.spct)
}


# untag -------------------------------------------------------------------


#' Remove tags
#'
#' Remove tags from an R object if present, otherwise return the object
#' unchanged.
#'
#' @param x an R object
#' @param ... not used in current version
#' @export untag
#'
#' @family tagging and related functions
#'
untag <- function(x, ...) UseMethod("untag")

#' @describeIn untag Default for generic function
#'
#' @export
#'
untag.default <- function(x, ...) {
  return(x)
}

#' @describeIn untag Specialization for generic_spct
#'
#' @param byref logical indicating if new object will be created by reference or
#'   by copy of x
#'
#' @return if \code{x} contains tag data they are removed and the "spct.tags"
#'   atrribute is set to \code{NA}, while if \code{x} has no tags, it is not
#'   modified. In either case, the byref argument is respected: in all cases if
#'   \code{byref=FALSE} a copy of \code{x} is returned.
#'
#' @export
#'
untag.generic_spct <- function(x,
                               byref = FALSE, ...) {
  if (!byref) {
    x <- x
    name <- NA
  } else {
    name <- substitute(x)
  }
  if (!is_tagged(x)) {
    return(x)
  }
  x[["wl.color"]] <- NULL
  x[["wb.f"]] <- NULL
  tag.data <- NA
  attr(x, "spct.tags") <- tag.data
  # to work by reference we need to assign the new spct to the old one
  if (byref & is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}
