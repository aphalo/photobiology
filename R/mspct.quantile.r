#' Quantiles of a collection of spectra
#'
#' Method to compute the "parallel" quantiles of values across members of a
#' collection of spectra or of a spectral object containing multiple spectra in
#' long form.
#'
#' @details Method specializations compute the qunatiles at each wavelength
#'   across a group of spectra stored in an object of one of the classes defined
#'   in package 'photobiology'. Omission of NAs is done separately at each
#'   wavelength. Interpolation is not applied, so all spectra in \code{x} must
#'   share the same set of wavelengths. An error is triggered if this condition
#'   is not fulfilled.
#'
#' @param x An R object. Currently this package defines methods for collections
#'   of spectral objects.
#' @param probs numeric vector of probabilities with values in \eqn{[0, 1]}.
#' @param na.rm	logical. A value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param ...	Further arguments passed by name to \code{\link[stats]{quantile}}.
#' @param simplify logical If \code{TRUE} and a single quantile is computed,
#'   return an spectrum by itself instead of as a single member of a collection.
#'
#' @return If \code{x} is a collection spectral of objects, such as a
#'   \code{"filter_mspct"} object, the returned object is of same class as the
#'   the collection, such as \code{"filter_mspct"}, containing one member
#'   summary spectrum for each value in \code{probs}, with variables with names
#'   tagged for summaries other than mean, median or quantiles. If a single
#'   quantile is computed and \code{simplify = TRUE} a single spectrum such as
#'    \code{"filter_spct"} is returned.
#'
#' @note Objects of classes \code{raw_spct} and \code{cps_spct} can contain data
#'   from multiple scans in multiple variables or "columns". The methods accept
#'   as arguments objects of these classes only if spectra contain data for a
#'   single spectrometer scan. In the case of \code{cps_spct} objects, a single
#'   column can also contain data from multiple scans spliced into a single
#'   variable.
#'
#' @seealso See \code{\link[stats]{quantile}} for the \code{qunatile()} method
#'   used for the computations. Additional arguments can be paased by name to
#'   be forwarded to \code{quantile}.
#'
#' @export
#'
#' @examples
#' s_quantile(sun_evening.mspct)
#'
s_quantile <- function(x, probs, na.rm, ...) UseMethod("s_quantile")


#' @describeIn s_quantile
#'
#' @export
#'
s_quantile.default <-
  function(x, probs = NA, na.rm = FALSE, ...) {
    warning("Method 's_quantile()' not implementd for objects of class ",
            class(x)[1], ".")
    ifelse(is.any_mspct(x), do.call(class(x[[1]])[1], args = list()), NA)
  }

#' @describeIn s_quantile
#'
#' @export
#'
s_quantile.generic_spct <-
  function(x, probs = c(0.25, 0.5, 0.75), na.rm = FALSE, ..., simplify = TRUE) {
    if (getMultipleWl(x) > 1) {
      s_quantile(subset2mspct(x),
                 probs = probs, na.rm = na.rm, simplify = simplify, ...)
    } else {
      x
    }
  }

#' @describeIn s_quantile
#'
#' @export
#'
s_quantile.source_mspct <-
  function(x, probs = c(0.25, 0.5, 0.75), na.rm = FALSE, ..., simplify = TRUE) {
    z <- source_mspct()
    for (p in probs) {
      q.name <- paste("p", p, sep = "_")
      z[[q.name]] <-
        rowwise_source(x = x,
                       .fun = stats::quantile,
                       probs = p, na.rm = na.rm, names = FALSE,
                       .fun.name = paste("Quantile P=", p, " of", sep = ""))
    }
    if (simplify && length(z) == 1L) {
      z[[1]]
    } else {
      z
    }
}

#' @describeIn s_quantile
#'
#' @export
#'
s_quantile.response_mspct <-
  function(x, probs = c(0.25, 0.5, 0.75), na.rm = FALSE, ..., simplify = TRUE) {
    z <- response_mspct()
    for (p in probs) {
      q.name <- paste("p", p, sep = "_")
      z[[q.name]] <-
        rowwise_response(x = x,
                         .fun = stats::quantile,
                         probs = p, na.rm = na.rm,
                         .fun.name = paste("Quantile P=", p, " of", sep = ""))
    }
    if (simplify && length(z) == 1L) {
      z[[1]]
    } else {
      z
    }
  }


#' @describeIn s_quantile
#'
#' @export
#'
s_quantile.filter_mspct <-
  function(x, probs = c(0.25, 0.5, 0.75), na.rm = FALSE, ..., simplify = TRUE) {
    z <- filter_mspct()
    for (p in probs) {
      q.name <- paste("p", p, sep = "_")
      z[[q.name]] <-
        rowwise_filter(x = x,
                       .fun = stats::quantile,
                       probs = p, na.rm = na.rm,
                       .fun.name = paste("Quantile P=", p, " of", sep = ""))
    }
    if (simplify && length(z) == 1L) {
      z[[1]]
    } else {
      z
    }
  }

#' @describeIn s_quantile
#'
#' @export
#'
s_quantile.reflector_mspct <-
  function(x, probs = c(0.25, 0.5, 0.75), na.rm = FALSE, ..., simplify = TRUE) {
    z <- reflector_mspct()
    for (p in probs) {
      q.name <- paste("p", p, sep = "_")
      z[[q.name]] <-
        rowwise_reflector(x = x,
                          .fun = stats::quantile,
                          probs = p, na.rm = na.rm,
                          .fun.name = paste("Quantile P=", p, " of", sep = ""))
    }
    if (simplify && length(z) == 1L) {
      z[[1]]
    } else {
      z
    }
  }

#' @describeIn s_quantile
#'
#' @export
#'
s_quantile.calibration_mspct <-
  function(x, probs = c(0.25, 0.5, 0.75), na.rm = FALSE, ..., simplify = TRUE) {
    z <- calibration_mspct()
    for (p in probs) {
      q.name <- paste("p", p, sep = "_")
      z[[q.name]] <-
        rowwise_calibration(x = x,
                            .fun = stats::quantile,
                            probs = p, na.rm = na.rm,
                            .fun.name = paste("Quantile P=", p, " of", sep = ""))
    }
    if (simplify && length(z) == 1L) {
      z[[1]]
    } else {
      z
    }
  }

#' @describeIn s_quantile
#'
#' @export
#'
s_quantile.cps_mspct <-
  function(x, probs = c(0.25, 0.5, 0.75), na.rm = FALSE, ..., simplify = TRUE) {
    z <- cps_mspct()
    for (p in probs) {
      q.name <- paste("p", p, sep = "_")
      z[[q.name]] <-
        rowwise_cps(x = x,
                    .fun = stats::quantile,
                    probs = p, na.rm = na.rm,
                    .fun.name = paste("Quantile P=", p, " of", sep = ""))
    }
    if (simplify && length(z) == 1L) {
      z[[1]]
    } else {
      z
    }
  }

#' @describeIn s_quantile
#'
#' @export
#'
s_quantile.raw_mspct <-
  function(x, probs = c(0.25, 0.5, 0.75), na.rm = FALSE, ..., simplify = TRUE) {
    z <- raw_mspct()
    for (p in probs) {
      q.name <- paste("p", p, sep = "_")
      z[[q.name]] <-
        rowwise_raw(x = x,
                    .fun = stats::quantile,
                    probs = p, na.rm = na.rm,
                    .fun.name = paste("Quantile P=", p, " of", sep = ""))
    }
    if (simplify && length(z) == 1L) {
      z[[1]]
    } else {
      z
    }
  }
