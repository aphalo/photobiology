#' Trim (or expand) tails of the spectrum based on wavelength limits,
#' interpolating the values at the limits.
#'
#' Trimming is needed for example to remove short wavelength noise
#' when the measured spectrum extends beyond the known emission
#' spectrum of the measured light source. Occasionally one may
#' want also to expand the wavelength range.
#'
#' @usage trim_spct(spct, band=NULL, low.limit=min(spct), high.limit=max(spct), use.hinges=TRUE, fill=NULL)
#'
#' @param spct an object of class "generic.spct"
#' @param band a numeric vector of length two, or any other object for which function range() will return two
#' @param low.limit shortest wavelength to be kept (defaults to shortest w.length value)
#' @param high.limit longest wavelength to be kept (defaults to longest w.length value)
#' @param use.hinges logical, if TRUE (the default)
#' wavelengths in nm.
#' @param fill if fill==NULL then tails are deleted, otherwise tails or s.irrad are filled with the value of fill
#'
#' @return a spectrum of same class as input with its tails trimmed or expanded
#'
#' @note When expanding an spectrum, if fill==NULL, then expansion is not performed.
#' Band can be "wave_band" object, a numeric vector or a list of numeric vectors, or any other user-defined or built-in
#' object for which range() returns a numeric vector of legth two, that can be interpreted as wavelengths expressed in nm.
#'
#' @keywords manip misc
#' @export
#' @examples
#' head(sun.spct)
#' head(trim_spct(sun.spct, low.limit=300))
#' head(trim_spct(sun.spct, low.limit=300, fill=NULL))
#' head(trim_spct(sun.spct, low.limit=300, fill=NA))
#' head(trim_spct(sun.spct, low.limit=300, fill=0.0))
#' head(trim_spct(sun.spct, low.limit=100, fill=0.0))
#' tail(trim_spct(sun.spct, low.limit=300, high.limit=1000, fill=NA))
#' tail(trim_spct(sun.spct, low.limit=300, high.limit=1000, fill=0.0))
#' tail(trim_spct(sun.spct, low.limit=300, high.limit=1000))
#' head(trim_spct(sun.spct, low.limit=300, high.limit=1000))
#' tail(trim_spct(sun.spct, low.limit=300, high.limit=400, fill=NA))
#' tail(trim_spct(sun.spct, low.limit=100, high.limit=400, fill=0.0))
#' head(trim_spct(sun.spct, low.limit=100, high.limit=400, fill=0.0))
#' head(trim_spct(sun.spct, band=new_waveband(300, 350)))
#' head(trim_spct(sun.spct, band=c(300, 350)))
#'

trim_spct <- function(spct, band=NULL, low.limit=min(spct), high.limit=max(spct), use.hinges=TRUE, fill=NULL)
{
  verbose <- FALSE # use option in the future

  names.spct <- names(spct)
  names.data <-names.spct[names.spct != "w.length"]
  class.spct  <- class(spct)
  comment.spct <- comment(spct)
  # check for target
  if (!is.null(band)) {
    trim.range <- range(band)
    low.limit <- trim.range[1]
    high.limit <- trim.range[2]
  }
  # check whether we should expand the low end
  low.end <- min(spct)
  if (low.end > low.limit) {
    if (!is.null(fill)) {
      # expand short tail
      low.tail.length <- low.end - low.limit
      low.tail.w.length <- seq(from = low.limit, to = low.end - 1, length=low.tail.length)
      spct.top <- data.table(w.length = low.tail.w.length)
      setattr(spct.top, "class", class.spct)
      for (data.col in names.data) {
        spct.top[ , eval(data.col) := fill]
      }
      spct <- rbind(spct.top, spct)
      setattr(spct, "comment", comment.spct)
      low.end <- min(spct)
    } else {
      if (verbose) {
        warning("Not trimming short end as low.limit is outside spectral data range.")
      }
      low.limit <- low.end
    }
  }

  # check whether we should expand the high end
  high.end <- max(spct)
  if (high.end < high.limit) {
    if (!is.null(fill)) {
      # expand short tail
      high.tail.length <- high.limit - high.end
      high.tail.w.length <- seq(from = high.end + 1, to = high.limit, length = high.tail.length)
      spct.bottom <- data.table(w.length = high.tail.w.length)
      setattr(spct.bottom, "class", class.spct)
      for (data.col in names.data) {
        spct.bottom[ , eval(data.col) := fill]
      }
      spct <- rbind(spct, spct.bottom)
      setattr(spct, "comment", comment.spct)
      low.end <- max(spct)
    } else {
      if (verbose) {
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
    return(spct[w.length %between% trim.range])
  }
  else {
    for (data.col in names.data) {
      spct[!w.length %between% trim.range, eval(data.col) := fill]
    }
  return(spct)
  }
}
