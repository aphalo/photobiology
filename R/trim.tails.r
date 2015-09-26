#' Trim (or expand) head and/or tail
#'
#' Trimming of tails of a spectrum based on wavelength limits, interpolating the
#' values at the boundaries.Trimming is needed for example to remove short
#' wavelength noise when the measured spectrum extends beyond the known emission
#' spectrum of the measured light source. Occasionally one may want also to
#' expand the wavelength range.
#'
#' @param w.length numeric array of wavelengths (nm)
#' @param s.irrad numeric array of spectral irradiance values
#' @param low.limit shortest wavelength to be kept (defaults to shortest
#'   w.length value)
#' @param high.limit longets wavelength to be kept (defaults to longest w.length
#'   value)
#' @param use.hinges logical, if TRUE (the default)
#' @param fill if fill==NULL then tails are deleted, otherwise tails or s.irrad
#'   are filled with the value of fill
#'
#' @return a data.frame with variables \code{w.length} and \code{s.irrad}
#'
#' @note When expanding an spectrum, if fill==NULL, then expansion is not
#'   performed with a warning.
#' @keywords manip misc
#' @family trim functions
#' @export
#' @examples
#' head(sun.data)
#' head(with(sun.data,
#'      trim_tails(w.length, s.e.irrad, low.limit=300)))
#' head(with(sun.data,
#'      trim_tails(w.length, s.e.irrad, low.limit=300, fill=NULL)))
#'
trim_tails <- function(w.length, s.irrad,
                       low.limit=min(w.length), high.limit=max(w.length),
                       use.hinges=TRUE, fill=NULL)
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
      low.limit <- low.end
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
      high.limit <- high.end
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
    return(dplyr::data_frame(w.length = w.length[trimmed.selector], s.irrad = s.irrad[trimmed.selector]))
  }
  else {
    s.irrad[!trimmed.selector] <- fill
    return(dplyr::data_frame(w.length = w.length, s.irrad=s.irrad))
  }
}
