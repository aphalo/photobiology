#' Get peaks in a spectrum
#'
#' This function finds all peaks (local maxima) in a spectrum, using a user
#' selectable size threshold relative to the tallest peak (global maximum). This
#' a wrapper built on top of function peaks from package splus2R.
#'
#' @param x numeric
#' @param y numeric
#' @param ignore_threshold numeric Value between 0.0 and 1.0 indicating the
#'   relative size compared to talelst peakthreshold below which peaks will be
#'   ignored.
#' @param span numeric A peak is defined as an element in a sequence which is
#'   greater than all other elements within a window of width \code{span}
#'   centered at that element. For example, a value of 3 means that a peak is
#'   bigger than both of its neighbors.
#' @param strict logical Flag: if TRUE, an element must be strictly greater than
#'   all other values in its window to be considered a peak. Default: TRUE.
#' @param x_unit character Vector of texts to be pasted at end of labels built
#'   from x value at peaks.
#' @param x_digits numeric Number of significant digits in wavelength label.
#'
#' @return A dataframe with variables w.length and s.irrad with their values at
#'   the peaks plus a character variable of labels.
#' @keywords manip misc
#' @export
#' @examples
#' with(sun.spct, get_peaks(w.length, s.e.irrad))
#'
#' @family peaks and valleys functions
#'
get_peaks <- function(x, y,
                      ignore_threshold = 0.0,
                      span = 5,
                      strict = TRUE,
                      x_unit = "",
                      x_digits = 3) {
  selector <- find_peaks(y, ignore_threshold, span, strict)
  if (sum(selector) < 1) {
    return(data.frame(x=numeric(0), y=numeric(0), label=character(0)))
  }
  else {
    peaks.x <- x[selector]
    peaks.y <- y[selector]
    if (length(peaks.x) == length(peaks.y)) {
      return(data.frame(x=peaks.x, y=peaks.y,
                        label=paste(as.character(signif(x=peaks.x, digits=x_digits)), x_unit, sep="")))
    }
    else {
      return(NA)
    }
  }
}
