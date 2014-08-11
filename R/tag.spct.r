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
                             use.hinges=NULL,
                             short.names=TRUE,
                             byref=TRUE, ...) {
  if (!byref) {
    x <- copy(x)
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
  wb.number <- length(w.band) # number of wavebands in list
  wb.name <- names(w.band) # their names in the list
  if (is.null(wb.name)) {
    wb.name <- character(wb.number)
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
  wb.rgb <- character(wb.number)
  wb.wl.low <- wb.wl.high <- numeric(wb.number)
  i <- 0L
  for (wb in w.band) {
    i <- i + 1L
    if (wb.name[i] == "") {
      if (short.names) {
        wb.name[i] <- labels(wb)[["label"]]
      } else {
        wb.name[i] <- labels(wb)[["name"]]
      }
    }
    wb.wl.low[i] <- min(wb)
    wb.wl.high[i] <- max(wb)
    wb.rgb[i] <- color(wb)[1]
  }
  n <- i
  x[ , idx := as.integer(NA) ]
  for (i in 1L:n) {
    x[ w.length >= wb.wl.low[i] & w.length < wb.wl.high[i], idx := as.integer(i) ]
  }
  x[ , band.name := wb.name[idx] ]
  x[ , band.color := wb.rgb[idx] ]
  x[ band.name == "none", band.color := "#000000" ]
  x[ , band.f := factor(band.name, levels=wb.name) ]
  x[ , idx := NULL]
}
