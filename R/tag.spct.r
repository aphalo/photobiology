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
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @param short.names logical indicating whether to use short or long names for wavebands
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @param ... not used in current version
#' @export tag.generic.spct
#'
tag.generic.spct <- function(x,
                             w.band=NULL,
                             use.hinges=TRUE,
                             short.names=TRUE,
                             byref=TRUE, ...) {
  if (!byref) {
    x <- copy(x)
    name <- NA
  } else {
    name <- substitute(x)
  }
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
  # we add a waveband for the whole spectrum
#  w.band <- c(list(new_waveband(min(x), max(x) + 1e-4, wb.name=NA)), w.band)
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
  tag.data <- list(time.unit=attr(x, "time.unit"),
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
  invisible(return(x))
}

#' Specialization for source.spct
#'
#' Tag a source.spct object using a list of wavebands.
#'
#' @param x a source.spct object
#' @param w.band list of waveband definitions created with new_waveband()
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @param short.names logical indicating whether to use short or long names for wavebands
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @param ... not used in current version
#' @export tag.source.spct
#'
#' @note At the moment not doing anything
#'
tag.source.spct <- function(x,
                             w.band=NULL,
                             use.hinges=NULL,
                             short.names=TRUE,
                             byref=TRUE, ...) {
  if (!byref) {
    x <- copy(x)
  } else {
  name <- substitute(x)
  }
  tag.generic.spct(x, w.band, use.hinges, short.names, byref=TRUE)
  #  x[ , irrad.color := s_e_irrad2rgb(w.length, s.e.irrad)]
  #  setattr(x, "spct.tags", )
  if (byref && is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(return(x))
}
