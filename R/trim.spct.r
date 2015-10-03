#' Trim (or expand) head and/or tail
#'
#' Trimming of head and tail of a spectrum based on wavelength limits,
#' interpolating the values at the boundaries. Trimming is needed for example to
#' remove short wavelength noise when the measured spectrum extends beyond the
#' known emission spectrum of the measured light source. Occasionally one may
#' want also to expand the wavelength range.
#'
#' @param spct an object of class "generic_spct"
#' @param range a numeric vector of length two, or any other object for which
#'   function range() will return two
#' @param low.limit shortest wavelength to be kept (defaults to shortest
#'   w.length value)
#' @param high.limit longest wavelength to be kept (defaults to longest w.length
#'   value)
#' @param use.hinges logical, if TRUE (the default) wavelengths in nm.
#' @param fill if fill==NULL then tails are deleted, otherwise tails or s.irrad
#'   are filled with the value of fill
#' @param byref logical indicating if new object will be created by reference or
#'   by copy of spct
#' @param verbose logical
#'
#' @return a spectrum of same class as input with its tails trimmed or expanded
#'
#' @note When expanding an spectrum, if fill==NULL, then expansion is not
#'   performed. Range can be "waveband" object, a numeric vector or a list of
#'   numeric vectors, or any other user-defined or built-in object for which
#'   \code{range()} returns a numeric vector of legth two, that can be
#'   interpreted as wavelengths expressed in nm.
#' @family trim functions
#' @keywords manip misc
#' @export
#' @examples
#' trim_spct(sun.spct, low.limit=300)
#' trim_spct(sun.spct, low.limit=300)
#' trim_spct(sun.spct, low.limit=300, fill=NULL)
#' trim_spct(sun.spct, low.limit=300, fill=NA)
#' trim_spct(sun.spct, low.limit=300, fill=0.0)
#'
trim_spct <- function(spct, range=NULL, low.limit=NULL, high.limit=NULL,
                      use.hinges=TRUE, fill=NULL, byref=FALSE, verbose=TRUE)
{
  if (is.null(spct)) {
    return(spct)
  }
  if (is.numeric(range)) {
    stopifnot(length(range) > 1)
  }
  if (!is.null(range) &&
      (!is.numeric(range) || (is.numeric(range) && length(range) != 2))) {
    range <- range(range, na.rm = TRUE)
  }
  if (is.null(use.hinges)) {
    use.hinges <- auto_hinges(spct)
  }
  stopifnot(is.any_spct(spct))
  if (byref) {
    name <- substitute(spct)
  }
  class_spct <- class(spct)
  if (!is.null(range)) {
    if (length(range) == 2 && is.na(range[1])) {
      low.limit <- NULL
    } else {
      low.limit <- ifelse(!is.null(low.limit), max(min(range, na.rm = TRUE), low.limit),
                          min(range, na.rm = TRUE))
    }
    if (length(range) == 2 && is.na(range[2])) {
      high.limit <- NULL
    } else {
      high.limit <- ifelse(!is.null(high.limit), min(max(range, na.rm = TRUE), high.limit),
                           max(range, na.rm = TRUE))
    }
  }
  trim.low <- !is.null(low.limit)
  trim.high <- !is.null(high.limit)
  if (trim.low && trim.high && high.limit - low.limit < 1e-7) {
    warning("When trimming 'range' must be a finite wavelength interval > 1E-7 nm")
    return(spct[FALSE, ]) # returns a spct object with nrow equal to zero
  }
  names.spct <- names(spct)
  names.data <- names.spct[names.spct != "w.length"]
  comment.spct <- comment(spct)
  time.unit.spct <- getTimeUnit(spct)
  Tfr.type.spct <- getTfrType(spct)
  Rfr.type.spct <- getRfrType(spct)
  # check whether we should expand the low end
  low.end <- min(spct, na.rm = TRUE)
  if (trim.low && low.end > low.limit) {
    if (!is.null(fill)) {
      # expand short tail
      low.tail.length <-  trunc(low.end - low.limit) + ifelse(use.hinges, 2, 1)
      low.tail.w.length <- seq(from = low.limit,
                               to = ifelse(use.hinges, low.end - 1e-12, low.end - 1),
                               length = low.tail.length)
      spct.top <- dplyr::data_frame(w.length = low.tail.w.length)
      for (data.col in names.data) {
        spct.top[[data.col]] <- fill
      }
      spct <- plyr::rbind.fill(list(spct.top, spct))
      spct <- dplyr::as_data_frame(spct)
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
  if (trim.high && high.end < high.limit) {
    if (!is.null(fill)) {
      # expand short tail
      high.tail.length <- trunc(high.limit - high.end) + ifelse(use.hinges, 2, 1)
      high.tail.w.length <- seq(from = ifelse(use.hinges, high.end + 1e-12, high.end + 1),
                                to = high.limit,
                                length = high.tail.length)
      spct.bottom <- dplyr::data_frame(w.length = high.tail.w.length)
      for (data.col in names.data) {
        spct.bottom[[data.col]] <- fill
      }
      spct <- plyr::rbind.fill(list(spct, spct.bottom))
      spct <- dplyr::as_data_frame(spct)
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
  } else {

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
      spct[!within.selector, data.col] <- fill
    }
  }
  # we now use plyr::rbind.fill which does not remove attributes
  # most of the code below may be redundant!!!
  class(spct) <- class_spct
  if (!is.null(comment.spct)) {
    comment(spct) <- comment.spct
  }
  if (!is.null(time.unit.spct) && !is.na(time.unit.spct)) {
    setTimeUnit(spct, time.unit.spct)
  }
  if (!is.null(Tfr.type.spct)) {
    setTfrType(spct, Tfr.type.spct)
  }
  if (!is.null(Rfr.type.spct)) {
    setRfrType(spct, Rfr.type.spct)
  }
  if (byref && is.name(name)) {
    name <- as.character(name)
    assign(name, spct, parent.frame(), inherits = TRUE)
  }
  check(spct)
  spct
}

#' @rdname trim_spct
#'
#' @param mspct an object of class "generic_mspct"
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
                       verbose = TRUE) {
  name <- substitute(mspct)

  z <- msmsply(mspct = mspct,
               .fun = trim_spct,
               range = range,
               low.limit = low.limit,
               high.limit = high.limit,
               use.hinges = use.hinges,
               fill = fill,
               byref = FALSE,
               verbose = verbose )

  if (byref & is.name(name)) {
    name <- as.character(name)
    assign(name, z, parent.frame(), inherits = TRUE)
  }
  z
}
