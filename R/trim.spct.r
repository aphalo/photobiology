#' Trim (or expand) tails of the spectrum based on wavelength limits,
#' interpolating the values at the limits.
#'
#' Trimming is needed for example to remove short wavelength noise
#' when the measured spectrum extends beyond the known emission
#' spectrum of the measured light source. Occasionally one may
#' want also to expand the wavelength range.
#'
#' @usage trim_spct(spct, range=NULL, low.limit=NULL, high.limit=NULL,
#'                  use.hinges=TRUE, fill=NULL, byref=FALSE, verbose=TRUE)
#'
#' @param spct an object of class "generic.spct"
#' @param range a numeric vector of length two, or any other object for which function range() will return two
#' @param low.limit shortest wavelength to be kept (defaults to shortest w.length value)
#' @param high.limit longest wavelength to be kept (defaults to longest w.length value)
#' @param use.hinges logical, if TRUE (the default)
#' wavelengths in nm.
#' @param fill if fill==NULL then tails are deleted, otherwise tails or s.irrad are filled with the value of fill
#' @param byref logical indicating if new object will be created by reference or by copy of spct
#' @param verbose logical
#'
#' @return a spectrum of same class as input with its tails trimmed or expanded
#'
#' @note When expanding an spectrum, if fill==NULL, then expansion is not performed.
#' Range can be "waveband" object, a numeric vector or a list of numeric vectors, or any other user-defined or built-in
#' object for which range() returns a numeric vector of legth two, that can be interpreted as wavelengths expressed in nm.
#'
#' @keywords manip misc
#' @export
#' @examples
#' trim_spct(sun.spct, low.limit=300)
#' my.sun.spct <- copy(sun.spct)
#' trim_spct(my.sun.spct, low.limit=300, byref=TRUE)
#' my.sun.spct
#' trim_spct(sun.spct, low.limit=300, fill=NULL)
#' trim_spct(sun.spct, low.limit=300, fill=NA)
#' trim_spct(sun.spct, low.limit=300, fill=0.0)
#' trim_spct(sun.spct, low.limit=300, high.limit=1000, fill=NA)
#' trim_spct(sun.spct, low.limit=300, high.limit=1000, fill=0.0)
#' trim_spct(sun.spct, low.limit=300, high.limit=1000)
#' trim_spct(sun.spct, low.limit=300, high.limit=400, fill=NA)
#' trim_spct(sun.spct, low.limit=100, high.limit=400, fill=0.0)
#' trim_spct(sun.spct, range=new_waveband(300, 350))
#' trim_spct(sun.spct, range=c(300, 350))
#' trim_spct(sun.spct, range=c(NA, 350))
#' trim_spct(sun.spct, range=c(NA, 300, 350))
#'

trim_spct <- function(spct, range=NULL, low.limit=NULL, high.limit=NULL, use.hinges=TRUE, fill=NULL, byref=FALSE, verbose=TRUE)
{
  if (is.null(spct)) {
    return(spct)
  }
  stopifnot(is.any.spct(spct))
  if (byref) {
    name <- substitute(spct)
  }
  class.spct <- class(spct)
  if (!is.null(range)) {
    if (length(range) == 2 && is.na(range[1])) {
      low.limit <- NULL
    } else {
      low.limit <- ifelse(!is.null(low.limit), max(min(range, nna.rm = TRUE), low.limit), min(range, na.rm = TRUE))
    }
    if (length(range) == 2 && is.na(range[2])) {
      high.limit <- NULL
    } else {
      high.limit <- ifelse(!is.null(high.limit), min(max(range, na.rm = TRUE), high.limit), max(range, na.rm = TRUE))
    }
  }
  if (is.null(low.limit)) {
    low.limit <- min(spct, na.rm=TRUE)
  }
  if (is.null(high.limit)) {
    high.limit <- max(spct, na.rm=TRUE)
  }
  if (high.limit - low.limit < 1e-7) {
    warning("When trimming 'range' must be a finite wavelength interval")
    return(NA) # this should be replaced with an empty spct object
  }
  names.spct <- names(spct)
  names.data <- names.spct[names.spct != "w.length"]
  comment.spct <- comment(spct)
  time.unit.spct <- getTimeUnit(spct)
  Tfr.type.spct <- getTfrType(spct)
  # check whether we should expand the low end
  low.end <- min(spct, na.rm=TRUE)
  if (low.end > low.limit) {
    if (!is.null(fill)) {
      # expand short tail
      low.tail.length <- low.end - low.limit
      low.tail.w.length <- seq(from = low.limit, to = low.end - 1, length=low.tail.length)
      spct.top <- data.table(w.length = low.tail.w.length)
      for (data.col in names.data) {
        spct.top[ , eval(data.col) := fill]
      }
      spct <- rbindlist(list(spct.top, spct))
      setGenericSpct(spct)
      low.end <- min(spct)
    } else {
      if (verbose) {
        # give a warning only if difference is > 0.01 nm
        if (verbose && (low.end - low.limit) > 0.01) {
          warning("Not trimming short end as low.limit is outside spectral data range.")
        }
      }
      low.limit <- low.end
    }
  }

  # check whether we should expand the high end
  high.end <- max(spct, na.rm=TRUE)
  if (high.end < high.limit) {
    if (!is.null(fill)) {
      # expand short tail
      high.tail.length <- high.limit - high.end
      high.tail.w.length <- seq(from = high.end + 1, to = high.limit, length = high.tail.length)
      spct.bottom <- data.table(w.length = high.tail.w.length)
      for (data.col in names.data) {
        spct.bottom[ , eval(data.col) := fill]
      }
      spct <- rbindlist(list(spct, spct.bottom))
      setGenericSpct(spct)
      low.end <- max(spct)
    } else {
      # give a warning only if difference is > 0.01 nm
      if (verbose && (high.limit - high.end) > 0.01) {
        warning("Not trimming long end as high.limit is outside spectral data range.")
      }
      high.limit <- high.end
    }
  }
  trim.range <- c(low.limit, high.limit)

  # insert hinges
  if (use.hinges) {
    hinges <- c(low.limit - 1e-4, low.limit, high.limit, high.limit + 1e-4)
    spct <- insert_spct_hinges(spct, hinges)
  }
  setkey(spct, w.length)
  if (is.null(fill)) {
    spct <- spct[w.length %between% trim.range]
  }
  else {
    for (data.col in names.data) {
      spct[!w.length %between% trim.range, eval(data.col) := fill]
    }
  }
  # we use rbindlist which removes derived class attributes
  setattr(spct, "class", class.spct)
  if (!is.null(comment.spct)) {
    setattr(spct, "comment", comment.spct)
  }
  if (!is.null(time.unit.spct)) {
    setTimeUnit(spct, time.unit.spct)
  }
  if (!is.null(Tfr.type.spct)) {
    setTfrType(spct, Tfr.type.spct)
  }
  if (byref && is.name(name)) {
    name <- as.character(name)
    assign(name, spct, parent.frame(), inherits = TRUE)
  }
  check(spct)
  return(spct)
}
