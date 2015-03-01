#' Generic function
#'
#' Tag values in an R object contains the expected data members.
#'
#' @param x an R object
#' @param ... not used in current version
#' @export tag
tag <- function(x, ...) UseMethod("tag")

#' Default for generic function
#'
#' Tag values in an R object contains the expected data members.
#'
#' @param x an R object
#' @param ... not used in current version
#' @export tag.default
tag.default <- function(x, ...) {
  return(x)
}

#' Specialization for generic.spct
#'
#' Tag a generic.spct object using a list of wavebands.
#'
#' @param x a generic.spct object
#' @param w.band list of waveband definitions created with new_waveband()
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @param short.names logical indicating whether to use short or long names for wavebands
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @param ... not used in current version
#' @export tag.generic.spct
#'
tag.generic.spct <- function(x,
                             w.band=NULL,
                             wb.trim=NULL,
                             use.hinges=TRUE,
                             short.names=TRUE,
                             byref=TRUE, ...) {
  if (!byref) {
    x <- copy(x)
    name <- NA
  } else {
    name <- substitute(x)
  }
  if (is.tagged(x)) {
    warning("Overwriting old tags in spectrum")
    untag(x)
  }
#   # we add a waveband for the whole spectrum
#   if (is.null(w.band)) {
#     w.band <- waveband(range(x))
#   }
  if (!is.null(w.band) && is.na(w.band[1])) {
    x[ , wl.color := w_length2rgb(x$w.length)]
    tag.data <- list(wl.color=TRUE)
    setattr(x, "spct.tags", tag.data)
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
  w.band <- trim_waveband(w.band=w.band, range=x, trim=wb.trim)
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
    length.wl <- length(x$w.length)
    use.hinges <- (x$w.length[length.wl] - x$w.length[1]) / length.wl > 0.2
    # we use 1.1 nm as performance degradation by using hinges is very significant
    # in the current version.
  }
  # we collect all hinges and insert them in one go
  # this may alter a little the returned values
  # but should be faster
  if (use.hinges) {
    all.hinges <- NULL
    for (wb in w.band) {
      if (!is.null(wb$hinges) & length(wb$hinges)>0) {
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
        wbs.name[i] <- ifelse(grepl("^range.", name.temp, ignore.case=TRUE),
                              paste("wb", i, sep=""),
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
  x[ , idx := n + 1L ]
  for (i in 1L:n) {
    x[ w.length >= wbs.wl.low[i] & w.length < wbs.wl.high[i], idx := as.integer(i) ]
  }
  wl.color.tmp <- w_length2rgb(x$w.length)
  x[ , wl.color := wl.color.tmp]
  x[ , wb.f := factor(wbs.name[idx], levels=wbs.name) ]
  x[ , idx := NULL]
  tag.data <- list(time.unit=getTimeUnit(x),
                   wb.key.name="Bands",
                   wl.color=TRUE,
                   wb.color=TRUE,
                   wb.num = n,
                   wb.colors=wbs.rgb[1:n],
                   wb.names=wbs.name[1:n],
                   wb.list=w.band)
  setattr(x, "spct.tags", tag.data)
  # to work by reference we need to assign the new DT to the old one
  if (byref & is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}


#' Make spectrum from a list of wavebands
#'
#' Make a generic.spct object with wavelengths from the range of wavebands
#' in a list.
#'
#' @param w.band list of waveband definitions created with new_waveband() or a single waveband object
#' @export
#'
#' @return a generic.spectrum object, with columns w.length as s.e.irrad,
#' the second one, se to 1 for all wavelengths.
#'
wb2spct <- function(w.band) {
  if (is(w.band, "waveband")) {
    w.band <- list(w.band)
  }
  w.length <- numeric(0)
  for (wb in w.band) {
    if (is(wb, "waveband")) {
      w.length <- c(w.length, wb$hinges)
    }
  }
  if (is.null(w.length) || length(w.length) < 2) {
    return(NA)
  }
  w.length <- unique(sort(w.length))
  new.spct <- data.table(w.length = w.length, s.e.irrad = 0, s.q.irrad = 0, Tfr = 0, Rfl = 0, s.e.response = 0)
  setGenSpct(new.spct)
  return(new.spct)
}

#' Make a tagged generic spectrum from a list of wavebands
#'
#' Make a tagged generic.spct object with wavelengths from the range of wavebands
#' in a list, and names of the same bands as factor levels, and also corresponding
#' colours.
#'
#' @param w.band list of waveband definitions created with new_waveband() or a single waveband object
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @param short.names logical indicating whether to use short or long names for wavebands
#' @param ... not used in current version
#' @export
#'
#' @note At the moment not doing anything
#'
#' @return a tagged spectrum
#'
wb2tagged_spct <- function(w.band,
                    use.hinges=TRUE,
                    short.names=TRUE,
                    ...) {
  new.spct <- wb2spct(w.band)
  tag(new.spct, w.band, use.hinges, short.names, byref=TRUE)
  new.spct[ , y := 0]
  return(new.spct)
}

#' Make spectrum from a list of wavebands
#'
#' Make a generic.spct object with wavelengths from the range of wavebands
#' in a list.
#'
#' @param w.band list of waveband definitions created with new_waveband() or a single waveband object
#' @param short.names logical indicating whether to use short or long names for wavebands
#' @export
#'
#' @return a generic.spectrum object, with columns w.length as s.e.irrad,
#' the second one, se to 1 for all wavelengths.
#'
wb2rect_spct <- function(w.band,
                         short.names=TRUE) {
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
        wbs.name[i] <- ifelse(grepl("^range.", name.temp, ignore.case=TRUE),
                              paste("wb", i, sep=""),
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
  new.spct <- data.table(w.length = wbs.wl.mid,
                         s.e.irrad = 0, s.q.irrad = 0, Tfr = 0, Rfl = 0, s.e.response = 0,
                         wl.color = w_length2rgb(wbs.wl.mid),
                         wb.f = factor(wbs.name, levels=wbs.name),
                         wl.high = wbs.wl.high, wl.low = wbs.wl.low,
                         y = 0)
  setGenSpct(new.spct)
  tag.data <- list(time.unit="none",
                   wb.key.name="Bands",
                   wl.color=TRUE,
                   wb.color=TRUE,
                   wb.num = wbs.number,
                   wb.colors=wbs.rgb,
                   wb.names=wbs.name,
                   wb.list=w.band)
  setattr(new.spct, "spct.tags", tag.data)

  return(new.spct)
}
