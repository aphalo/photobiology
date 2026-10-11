#' Normalize a range argument into a true numeric range
#'
#' Several functions in this package and the suite accept a range argument
#' with a flexible syntax. To ensure that all functions and methods behave
#' in the same way this code has been factored out into a separate function.
#'
#' @param arg.range a numeric vector of length two, or any other object for
#'   which function range() will return a range of wavelengths (nm).
#' @param wl.range a numeric vector of length two, or any other object for which
#'   function range() will return a range of wavelengths (nm), missing values
#'   are not allowed.
#' @param trim logical If TRUE the range returned is bound within
#'   \code{wl.range} while if FALSE it can be broader.
#'
#' @return a numeric vector of length two, guaranteed not to have missing
#'   values.
#'
#' @details The \code{arg.range} argument can contain NAs which are replaced by
#'   the value at the same position in \code{wl.range}. In addition
#'   a NULL argument for \code{range} is converted into \code{wl.range}.
#'   The \code{wl.range} is also the limit to which the returned value
#'   is trimmed if \code{trim == TRUE}. The idea is that the value supplied as
#'   \code{wl.range} is the wavelength range of the data.
#'
#' @family auxiliary functions
#'
#' @export
#' @examples
#' normalize_range_arg(c(NA, 500), range(sun.spct))
#' normalize_range_arg(c(300, NA), range(sun.spct))
#' normalize_range_arg(c(100, 5000), range(sun.spct), FALSE)
#' normalize_range_arg(c(NA, NA), range(sun.spct))
#' normalize_range_arg(c(NA, NA), sun.spct)
#'
normalize_range_arg <- function(arg.range, wl.range, trim = TRUE) {
  if (!is.numeric(wl.range) ||
        (is.numeric(wl.range) && length(wl.range) != 2)) {
    wl.range <- range(wl.range)
  }
  stopifnot(is.numeric(wl.range) && length(unique(wl.range)) == 2)

  if (is.null(arg.range) || all(is.na(arg.range))) {
    return(wl.range)
  }
  if (!is.numeric(arg.range) ||
        (is.numeric(arg.range) && length(arg.range) != 2)) {
    arg.range <- range(arg.range, na.rm = TRUE)
  }
  stopifnot(is.numeric(arg.range) && length(arg.range) == 2)

  if (is.na(arg.range[1]) || trim && arg.range[1] < wl.range[1])
    arg.range[1] <- wl.range[1]

  if (is.na(arg.range[2]) || trim && arg.range[2] > wl.range[2])
    arg.range[2] <- wl.range[2]

  # NAs have been replaced above
  if (diff(arg.range) < 1e-3) {
    c(1, 2) # nm
  } else {
    arg.range
  }
}
