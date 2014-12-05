#' Trim (or expand) tails of the spectrum based on wavelength limits,
#' interpolating the values at the limits.
#'
#' Trimming is needed for example to remove short wavelength noise
#' when the measured spectrum extends beyond the known emission
#' spectrum of the measured light source. Occasionally one may
#' want also to expand the wavelength range.
#'
#' @usage trim_spct(spct, band=NULL, low.limit=NULL, high.limit=NULL, use.hinges=TRUE, fill=NULL, byref=FALSE)
#'
#' @param spct an object of class "generic.spct"
#' @param band a numeric vector of length two, or any other object for which function range() will return two
#' @param low.limit shortest wavelength to be kept (defaults to shortest w.length value)
#' @param high.limit longest wavelength to be kept (defaults to longest w.length value)
#' @param use.hinges logical, if TRUE (the default)
#' wavelengths in nm.
#' @param fill if fill==NULL then tails are deleted, otherwise tails or s.irrad are filled with the value of fill
#' @param byref logical indicating if new object will be created by reference or by copy of spct
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
#' trim_spct(sun.spct, band=new_waveband(300, 350))
#' trim_spct(sun.spct, band=c(300, 350))
#'

trim_spct <- function(spct, band=NULL, low.limit=NULL, high.limit=NULL, use.hinges=TRUE, fill=NULL, byref=FALSE)
{
  if (!is.any.spct(spct)) {
    setGenericSpct(spct)
  }
  verbose <- TRUE
  if (byref) {
    name <- substitute(spct)
  }
  if (is.null(low.limit)) {
    low.limit <- min(spct, na.rm=TRUE)
  }
  if (is.null(high.limit)) {
    high.limit <- max(spct, na.rm=TRUE)
  }
  names.spct <- names(spct)
  names.data <- names.spct[names.spct != "w.length"]
  #  class.spct  <- class(spct)
  comment.spct <- comment(spct)
  time.unit.spct <- attr(spct, "time.unit", exact=TRUE)
  Tfr.type.spct <- attr(spct, "Tfr.type", exact=TRUE)
  # check for target
  if (!is.null(band)) {
    trim.range <- range(band)
    low.limit <- trim.range[1]
    high.limit <- trim.range[2]
  }
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
      setattr(spct.top, "class", class(spct))
      spct <- rbindspct(list(spct.top, spct))
      low.end <- min(spct)
    } else {
      if (verbose) {
        # give a warning only if difference is > 0.01 nm
        if ((low.end - low.limit) > 0.01) {
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
      setattr(spct.bottom, "class", class(spct))
      spct <- rbindspct(list(spct, spct.bottom))
      low.end <- max(spct)
    } else {
      # give a warning only if difference is > 0.01 nm
      if ((high.limit - high.end) > 0.01) {
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
  if (!is.null(comment.spct)) {
    setattr(spct, "comment", comment.spct)
  }
  if (!is.null(time.unit.spct)) {
    setattr(spct, "time.unit", time.unit.spct)
  }
  if (!is.null(Tfr.type.spct)) {
    setattr(spct, "Tfr.typr", Tfr.type.spct)
  }
  if (byref && is.name(name)) {
    name <- as.character(name)
    assign(name, spct, parent.frame(), inherits = TRUE)
  }
  # %between% removes derived class tags!
  #  setattr(spct, "class", class.spct)
  invisible(spct)
}
