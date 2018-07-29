#' Fit peaks in a spectrum
#'
#' This function fits peaks (local maxima) in a spectrum.
#'
#' @param x, y numeric vectors
#' @param target.idxs integer Index into the vectors pointing to peak candidates.
#' @param span integer The number of observations to use in the fit of each
#'   peak.
#' @param method character.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before fittint the peaks.
#'
#' @return a list frame with three columns, \code{x, y, .
#'
#' @export
#'
#' @family peaks and valleys functions
#'
fit_peaks <-
  function(x,
           y,
           target.idxs = numeric(),
           span = 12,
           shape = "",
           na.rm = FALSE) {
    data <- tibble(x, y)
    if (na.rm) {
      data <- na.omit(data)
    }
    .fun.peak <- function(x, y, x.peak, y.peak, a, b, c, d) {
      xi <- x - x.peak
      y.peak + a * xi - b * ((xi + d)^c + 1)^(1/c)
    }
    for (i in target.idxs) {
      peak.data <- data[(i - span %/% 2L):(i + span %/% 2L), ]
      fm <- nls()
    }
    range_x <- range(x, finite = TRUE)
    min_x <- range_x[1]
    max_x <- range_x[2]
    x <- ifelse(!is.finite(x), min_x, x)
    # the next two lines cater for the case when max_x < 0, which is quite common with logs
    delta <- max_x - min_x
    top_flag <- ignore_threshold > 0.0
    scaled_threshold <- delta * abs(ignore_threshold)
    pks <- splus2R::peaks(x = x, span = span, strict = strict)
    if (abs(ignore_threshold) < 1e-5)
      return(pks)
    if (top_flag) {
      return(ifelse(x - min_x > scaled_threshold, pks , FALSE))
    } else {
      return(ifelse(max_x - x > scaled_threshold, pks , FALSE))
    }
  }

