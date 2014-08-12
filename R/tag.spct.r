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
  }
  if (!is.null(w.band) && is.na(w.band[1])) {
    x[ , wl.color := w_length2rgb(w.length)]
    return(x)
  }
  if (!is.null(w.band) & is(w.band, "waveband")) {
    # if the argument is a single w.band, we enclose it in a list
    # so that the for loop works as expected.This is a bit of a
    # cludge but lets us avoid treating it as a special case
    w.band <- list(w.band)
  }
  # we add a waveband for the whole spectrum
  w.band <- c(list(new_waveband(min(x), max(x) + 1e-4, wb.name="none")), w.band)
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
    use.hinges <- (x$w.length[length.wl] - x$w.length[1]) / length.wl > 1.1
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
        all.hinges <- c(all.hinges, wb$hinges)
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
        wbs.name[i] <- labels(wb)[["label"]]
      } else {
        wbs.name[i] <- labels(wb)[["name"]]
      }
    }
    wbs.wl.low[i] <- min(wb)
    wbs.wl.high[i] <- max(wb)
    wbs.rgb[i] <- color(wb)[1]
  }
  n <- i
  x[ , idx := as.integer(NA) ]
  for (i in 1L:n) {
    x[ w.length >= wbs.wl.low[i] & w.length < wbs.wl.high[i], idx := as.integer(i) ]
  }
  x[ , wl.color := w_length2rgb(w.length)]
  x[ , wb.color := wbs.rgb[idx] ]
  x[ , wb.name := wbs.name[idx] ]
  x[ wb.name == "none", wb.color := "#000000" ]
  x[ , wb.f := factor(wb.name, levels=wbs.name) ]
  x[ , idx := NULL]
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
tag.source.spct <- function(x,
                             w.band=NULL,
                             use.hinges=NULL,
                             short.names=TRUE,
                             byref=TRUE, ...) {
  if (!byref) {
    x <- copy(x)
  }
  q2e(x, byref=TRUE)
  x[ , irrad.color := s_e_irrad2rgb(w.length, s.e.irrad)]
  tag.generic.spct(x, w.band, use.hinges, short.names, byref)
}
