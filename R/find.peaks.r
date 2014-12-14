#' Find peaks in a spectrum.
#'
#' This function finds all peaks (local maxima) in a spectrum, using a user selectable 
#' size threshold relative to the tallest peak (global maximum). This a wrapper built
#' on top of function peaks from package splus2R. 
#'
#' @usage find_peaks(x, ignore_threshold=0.0, span=3, strict=TRUE)
#' @param x numeric array 
#' @param ignore_threshold numeric value between 0.0 and 1.0 indicating the size threshold below which peaks will be ignored.
#' @param span a peak is defined as an element in a sequence which is greater than all other elements within a window of width span centered at that element. The default value is 3, meaning that a peak is bigger than both of its neighbors. Default: 3.
#' @param strict ogical flag: if TRUE, an element must be strictly greater than all other values in its window to be considered a peak. Default: TRUE.
#' 
#' @return an object like s.irrad of logical values. Values that are TRUE correspond to local peaks in the data.
#' @keywords manip misc
#' @export
#' @examples
#' with(sun.data, w.length[find_peaks(s.e.irrad)])
#' @note this function handles non-finite (including NA) values differently than \code{peaks}, instead of giving an error
#' they are replaced with the smallest finite value in the \code{x}.
#' @importFrom splus2R peaks
#' 
find_peaks <- function(x, ignore_threshold=0.0, span=3, strict=TRUE) {
  range_x <- range(x, finite=TRUE)
  min_x <- range_x[1]
  max_x <- range_x[2]
  x <- ifelse(!is.finite(x), min_x, x)
  # the next two lines catter for the case when max_x < 0, which is quite common with logs
  delta <- max_x - min_x 
  top_flag <- ignore_threshold > 0.0
  scaled_threshold <- delta * abs(ignore_threshold)
  pks <- splus2R::peaks(x=x, span=span, strict=strict)
  if (abs(ignore_threshold) < 1e-5) return(pks)
  if (top_flag) {
    return(ifelse(x - min_x > scaled_threshold, pks , FALSE))
  } else {
    return(ifelse(max_x - x > scaled_threshold, pks , FALSE))
  }
}
