#' Trim tails of the spectrum based on wavelength limits,
#' interpolating the values at the limits.
#'
#' Trimming is needed for example to remove short wavelength noise
#' when the measured spectrum extends beyond the known emission
#' spectrum of the measured light source.
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
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' head(sun.data)
#' head(with(sun.data, trim_tails(w.length, s.e.irrad, low.limit=300)))

trim_tails <- function(w.length, s.irrad, low.limit=min(w.length), high.limit=max(w.length), use.hinges=TRUE, fill=NULL)
{
  # insert hinges
  if (use.hinges) {
    new.data <- insert_hinges(w.length, s.irrad, c(low.limit, high.limit))
    w.length <- new.data$w.length
    s.irrad <- new.data$s.irrad
  }
  trimmed.selector <- w.length >= low.limit & w.length <= high.limit
  if (is.null(fill)) {
    return(data.frame(w.length=w.length[trimmed.selector], s.irrad=s.irrad[trimmed.selector]))
  }
  else {
    return(data.frame(w.length=w.length, s.irrad=s.irrad[!trimmed.selector] <- fill))
  }
}
