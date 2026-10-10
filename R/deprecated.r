#' Get peaks and valleys from a spectrum
#'
#' These functions "get" (or extract) peaks (maxima) and valleys (minima) in two
#' vectors, usually a spectral quantity and wavelength, using a user selectable
#' span for window width and global and local (within moving window) size
#' thresholds. They also generate \code{character} values for \code{x}.
#'
#' @param x,y numeric
#' @param x_unit character Vector of texts to be pasted at end of labels built
#'   from x value at peaks.
#' @param x_digits numeric Number of significant digits in wavelength label.
#' @inheritParams find_peaks
#'
#' @inherit find_peaks details seealso
#'
#' @note The use of these two functions is deprecated. They are retained for
#'   backwards compatibility and will be removed in the near future.
#'
#' @return A data frame with variables w.length and s.irrad with their values at
#'   the peaks or valleys plus a character variable of labels.
#'
#' @export
#'
#' @family peaks and valleys functions
#'
get_peaks <- function(x,
                      y,
                      global.threshold = 0,
                      span = 5,
                      strict = TRUE,
                      x_unit = "",
                      x_digits = 3,
                      na.rm = FALSE) {
  message("Functions 'get_peaks()' and 'get_valeys()' have been deprecated: ",
          "please use 'peaks()' and 'valleys()', or 'find_peaks()' and ",
          "'find_valleys()', instead.")
  stopifnot(length(x) == length(y))
  selector <- find_peaks(x = y,
                         global.threshold = global.threshold,
                         span = span,
                         strict = strict,
                         na.rm = na.rm)
  if (sum(selector) < 1) {
    return(data.frame(
      x = numeric(0),
      y = numeric(0),
      label = character(0)
    ))
  } else {
    peaks.x <- x[selector]
    peaks.y <- y[selector]
    return(data.frame(
      x = peaks.x,
      y = peaks.y,
      label = paste(as.character(signif(
        x = peaks.x, digits = x_digits
      )), x_unit, sep = "")
    ))
  }
}

#' @rdname get_peaks
#' @export
#'
get_valleys <- function(x, y,
                        global.threshold = 0,
                        span = 5,
                        strict = TRUE,
                        x_unit = "",
                        x_digits = 3,
                        na.rm = FALSE) {
  xy.data <- get_peaks(x = x, y = -y,
                       global.threshold = global.threshold,
                       span = span,
                       strict = strict,
                       x_unit = x_unit,
                       x_digits = x_digits,
                       na.rm = na.rm)
  xy.data[["y"]] <- -xy.data[["y"]]
  return(xy.data)
}
