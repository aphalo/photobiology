#' Get peaks in a spectrum.
#'
#' This function finds all valleys (local minima) in a spectrum, using a user
#' selectable size threshold relative to the tallest peak (global maximum). This
#' a wrapper built on top of function peaks from package splus2R.
#'
#' @usage get_valleys(x, y, ignore_threshold = 0.0, span = 5, strict = TRUE,
#'   x_unit = "", x_digits = 3)
#' @param x numeric array
#' @param y numeric array
#' @param ignore_threshold numeric value between 0.0 and 1.0 indicating the
#'   relative size compared to talelst peakthreshold below which peaks will be
#'   ignored.
#' @param span a peak is defined as an element in a sequence which is greater
#'   than all other elements within a window of width span centered at that
#'   element. The default value is 3, meaning that a peak is bigger than both of
#'   its neighbors. Default: 3.
#' @param strict logical flag: if TRUE, an element must be strictly greater than
#'   all other values in its window to be considered a peak. Default: TRUE.
#' @param x_unit  character vector to be pasted at end of lables built from x
#'   value at peaks.
#' @param x_digits number of significant digits in wevelength label.
#'
#' @return a dataframe of variables w.length and s.irrad with their values at
#'   the peaks plus a character variable of labels.
#' @keywords manip misc
#' @export
#' @examples
#' with(sun.spct, get_valleys(w.length, s.e.irrad))
#'
#' @family peaks and valleys functions
#'
get_valleys <- function(x, y,
                        ignore_threshold = 0.0,
                        span = 5,
                        strict = TRUE,
                        x_unit = "",
                        x_digits = 3) {
  xy.data <- get_peaks(x, -y,
                       -ignore_threshold,
                       span = span,
                       strict = strict,
                       x_unit = x_unit,
                       x_digits = x_digits)
  xy.data$y <- -xy.data$y
  return(xy.data)
}
