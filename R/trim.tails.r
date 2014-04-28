#' Trim (or expand) tails of the spectrum based on wavelength limits,
#' interpolating the values at the limits.
#'
#' Trimming is needed for example to remove short wavelength noise
#' when the measured spectrum extends beyond the known emission
#' spectrum of the measured light source. Occasionally one may
#' want also to expand the wavelength range.
#' 
#' @usage trim_tails(w.length, s.irrad, low.limit=min(w.length), high.limit=max(w.length), use.hinges=TRUE, fill=NULL) 
#' 
#' @param w.length numeric array of wavelengths (nm)
#' @param s.irrad numeric array of spectral irradiance values
#' @param low.limit shortest wavelength to be kept (defaults to shortest w.length value)
#' @param high.limit longets wavelength to be kept (defaults to longest w.length value)
#' @param use.hinges logical, if TRUE (the default) 
#' @param fill if fill==NULL then tails are deleted, otherwise tails or s.irrad are filled with the value of fill
#' 
#' @return a data.frame with variables \code{w.length} and \code{s.irrad}
#' 
#' @note When expanding an spectrum, if fill==NULL, then expansion is not performed with a warning.
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' head(sun.data)
#' head(with(sun.data, trim_tails(w.length, s.e.irrad, low.limit=300)))
#' head(with(sun.data, trim_tails(w.length, s.e.irrad, low.limit=300, fill=NULL)))
#' head(with(sun.data, trim_tails(w.length, s.e.irrad, low.limit=300, fill=NA)))
#' head(with(sun.data, trim_tails(w.length, s.e.irrad, low.limit=300, fill=0.0)))
#' head(with(sun.data, trim_tails(w.length, s.e.irrad, low.limit=100, fill=0.0)))
#' \notrun{
#' tail(with(sun.data, trim_tails(w.length, s.e.irrad, low.limit=300, high.limit=1000, fill=NA)))
#' tail(with(sun.data, trim_tails(w.length, s.e.irrad, low.limit=300, high.limit=1000, fill=0.0)))
#' tail(with(sun.data, trim_tails(w.length, s.e.irrad, low.limit=300, high.limit=1000)))
#' head(with(sun.data, trim_tails(w.length, s.e.irrad, low.limit=300, high.limit=1000)))
#' tail(with(sun.data, trim_tails(w.length, s.e.irrad, low.limit=300, high.limit=400, fill=NA)))
#' tail(with(sun.data, trim_tails(w.length, s.e.irrad, low.limit=100, high.limit=400, fill=0.0)))
#' head(with(sun.data, trim_tails(w.length, s.e.irrad, low.limit=100, high.limit=400, fill=0.0)))
#' }

trim_tails <- function(w.length, s.irrad, low.limit=min(w.length), high.limit=max(w.length), use.hinges=TRUE, fill=NULL)
{
  if (!check_spectrum(w.length, s.irrad)) return(NA)

  # check whether we should expand the low end
  low.end <- min(w.length)
  if (low.end > low.limit) {
    if (!is.null(fill)) {
      # expand short tail
      low.tail.length <- low.end - low.limit
      low.tail <- seq(from = low.limit, to = low.end - 1, length=low.tail.length)
      s.irrad <- c(rep(fill, length.out = low.tail.length), s.irrad)
      w.length <- c(low.tail, w.length)
      low.end <- low.limit
    } else {
      warning("Ignoring low.limit as it is too low.")
      low.limit <- min(w.length)
    }
  }
  
  # check whether we should expand the high end
  high.end <- max(w.length)  
  if (high.end < high.limit) {
    if (!is.null(fill)) {
      # expand short tail
      high.tail.length <- high.limit - high.end
      high.tail <- seq(from = high.end + 1, to = high.limit, length = high.tail.length)
      s.irrad <- c(s.irrad, rep(fill, length.out = high.tail.length))
      w.length <- c(w.length, high.tail)
      high.end <- high.limit
    } else {
      warning("Ignoring high.limit as it is too high.")
      high.limit <- max(w.length)
    }
  }
  
  # insert hinges
  if (use.hinges) {
    new.data <- insert_hinges(w.length, s.irrad, c(low.limit, high.limit))
    w.length <- new.data$w.length
    s.irrad <- new.data$s.irrad
  }
  
  trimmed.selector <- (w.length >= low.limit) & (w.length <= high.limit)
  if (is.null(fill)) {
    return(data.frame(w.length = w.length[trimmed.selector], s.irrad = s.irrad[trimmed.selector]))
  }
  else {
    s.irrad[!trimmed.selector] <- fill
    return(data.frame(w.length = w.length, s.irrad=s.irrad))
  }
}
